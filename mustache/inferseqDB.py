import sys
import warnings
warnings.filterwarnings("ignore")
import click
import pandas as pd
from mustache import fastatools
from mustache import bowtie2tools
from mustache import sctools
from mustache import misc
from mustache import pysamtools
import pygogo as gogo
import pysam
from Bio import SeqIO
from collections import OrderedDict
from os.path import dirname, join
from random import randint
from snakemake import shell
from collections import defaultdict
from scipy.sparse.csgraph import connected_components
import numpy as np

verbose=True
logger = gogo.Gogo(__name__, verbose=verbose).logger


import warnings
warnings.filterwarnings("ignore", message="numpy.dtype size changed")
warnings.filterwarnings("ignore", message="numpy.ufunc size changed")


def _inferseq(pairsfile, bamfile, inferseq_assembly, inferseq_reference, output_file):

    index_genomes(inferseq_assembly, inferseq_reference)


    tmp_dir = dirname(output_file)

    pairs = pd.read_csv(pairsfile, sep='\t')

    logger.info("Aligning pairs to assembly...")
    assembly_flanks_fasta_prefix = write_flanks_to_align_to_assembly(pairs, bam, reference_genome_dict, tmp_dir)
    assembly_outbam = join(tmp_dir, 'mustache.inferseq.' + str(randint(0, 1e20)) + '.bam')
    bowtie2tools.align_paired_fasta_to_genome(
        assembly_flanks_fasta_prefix+'.R1.fasta',
        assembly_flanks_fasta_prefix+'.R2.fasta',
        inferseq_assembly, assembly_outbam, silence=True,
        additional_flags='-X 10000 --local -a --no-mixed --no-discordant --dovetail'
    )

    logger.info("Inferring sequences from pairs aligned to assembly and in context...")
    sequences_inferred_from_assembly_with_context = infer_sequences_with_context(assembly_outbam, assembly_genome_dict)
    logger.info("Inferring sequences from pairs aligned to assembly...")
    sequences_inferred_from_assembly_without_context = infer_sequences_without_context(assembly_outbam, assembly_genome_dict)

    shell('rm {fasta_prefix}* {outbam}*'.format(fasta_prefix=assembly_flanks_fasta_prefix, outbam=assembly_outbam))

    logger.info("Aligning pairs to reference genome...")
    reference_flanks_fasta_prefix = write_flanks_to_align_to_reference(pairs, tmp_dir)
    reference_outbam = join(tmp_dir, 'mustache.inferseq.' + str(randint(0, 1e20)) + '.bam')
    bowtie2tools.align_paired_fasta_to_genome(
        reference_flanks_fasta_prefix + '.R1.fasta',
        reference_flanks_fasta_prefix + '.R2.fasta',
        inferseq_reference, reference_outbam, silence=True,
        additional_flags='-X 10000 --local -a --no-mixed --no-discordant --dovetail'
    )

    logger.info("Inferring sequences from pairs aligned to reference...")
    sequences_inferred_from_reference = infer_sequences_without_context(
        reference_outbam, reference_genome_dict, softclip_filter='lenient', min_identity=0.95)

    shell('rm {fasta_prefix}* {outbam}*'.format(fasta_prefix=reference_flanks_fasta_prefix, outbam=reference_outbam))

    method1 = make_dataframe(sequences_inferred_from_assembly_with_context, method='inferred_assembly_with_context')
    method2 = make_dataframe(sequences_inferred_from_assembly_without_context, method='inferred_assembly_without_context')
    method3 = make_dataframe(sequences_inferred_from_reference, method='inferred_reference')

    all_inferred_results = method1.append(method2, ignore_index=True).append(method3, ignore_index=True).sort_values(
        by=['pair_id', 'method']
    )

    all_inferred_results.loc[:, 'pair_id'] = list(map(str, map(int, list(all_inferred_results['pair_id']))))
    all_inferred_results = all_inferred_results.query("inferred_seq_length > 0")

    logger.info("Writing results to file %s..." % output_file)
    all_inferred_results.to_csv(output_file, sep='\t', index=False)


def index_genomes(inferseq_assembly, inferseq_reference):
    logger.info("Indexing inferseq assembly and inferseq reference if not already indexed...")
    if not bowtie2tools.genome_is_indexed(inferseq_assembly):
        bowtie2tools.index_genome(inferseq_assembly)
    if not bowtie2tools.genome_is_indexed(inferseq_reference):
        bowtie2tools.index_genome(inferseq_reference)
    logger.info("Genome has been indexed...")
    
@click.command()
@click.argument('pairsfile', type=click.Path(exists=True))
@click.argument('fasta_database', type=click.Path(exists=True))
@click.option('--output_file', '-o', default='mustache.inferseq.tsv', help="The output file to save the results.")
def inferseq(pairsfile, fasta_database, output_file=None):
    _inferseq(pairsfile, fasta_database, output_file)



if __name__ == '__main__':
    inferseq()