import warnings
warnings.filterwarnings("ignore")
import sys
import click
import pysam
import pygogo as gogo
import pandas as pd
import numpy as np
from snakemake import shell
from random import randint
from mustache import fastatools, embosstools, sctools, pysamtools, minimustools
from mustache.misc import revcomp
from os.path import basename
from multiprocessing import Pool
from os.path import join, dirname

verbose=True
logger = gogo.Gogo(__name__, verbose=verbose).logger

def extend(row):
    bam = pysam.AlignmentFile(row['bam_path'], 'rb')
    contig, pos, orient, seq = row['contig'], row['pos'], row['orient'], row['consensus_seq']
    tmp_outdir = row['outdir']

    logger.info("Attempting to extend sequence at %s..." % ' '.join(map(str, [contig, pos, orient])))
    reads, quals = get_reads_to_assemble(bam, contig, pos, orient, get_quals=True)
    if len(reads) == 0:
        return seq

    #print("RUNNING ASSEMBLY")
    assembler = minimustools.MinimusAssembler(reads, quals, outdir=tmp_outdir)
    assembler.assemble()

    #print("CHECKING IF SOMETHING ASSEMBLED")
    if not assembler.something_assembled():
        assembler.delete_files()
        return seq

    #print("ALIGNING TO ASSEMBLY")
    assembler.align_seq_to_assembly(seq)
    #print("RETRIEVING EXTENDED SEQUENCE")
    extended_seq = assembler.retrieve_extended_sequence(orient)
    assembler.delete_files()

    if extended_seq is None:
        return seq

    logger.info("Sequence extended successfully by %d base pairs" % (len(extended_seq) - len(seq)))
    return extended_seq

def get_reads_to_assemble(bam, contig, pos, orient, get_quals=False):

    if get_quals:

        if orient == 'R':
            softclipped_reads, softclipped_quals = pysamtools.get_right_softclipped_reads_at_site(bam, contig, pos, get_quals=True)
            unmapped_reads, unmapped_quals = pysamtools.get_right_unmapped_reads(bam, contig, pos, get_quals=True)
        elif orient == 'L':
            softclipped_reads, softclipped_quals = pysamtools.get_left_softclipped_reads_at_site(bam, contig, pos, get_quals=True)
            unmapped_reads, unmapped_quals = pysamtools.get_left_unmapped_reads(bam, contig, pos, get_quals=True)

    else:
        if orient == 'R':
            softclipped_reads = pysamtools.get_right_softclipped_reads_at_site(bam, contig, pos, get_quals=True)
            unmapped_reads = pysamtools.get_right_unmapped_reads(bam, contig, pos, get_quals=True)
        elif orient == 'L':
            softclipped_reads = pysamtools.get_left_softclipped_reads_at_site(bam, contig, pos, get_quals=True)
            unmapped_reads = pysamtools.get_left_unmapped_reads(bam, contig, pos, get_quals=True)


    if get_quals:
        return softclipped_reads + unmapped_reads, softclipped_quals + unmapped_quals
    else:
        return softclipped_reads + unmapped_reads


def _extendflanks(flanksfile, bamfile, threads, output_file):
    flanks = pd.read_csv(flanksfile, sep='\t')

    if flanks.shape[0] == 0:
        logger.info("No flanks found in the input file...")

    else:
        sequences = list(flanks['consensus_seq'])
        did_extend = [False]*len(sequences)

        logger.info("Running extendflanks algorithm on %d total flanks..." % flanks.shape[0])

        flank_rows = [row for index, row in flanks.iterrows()]
        tmp_outdir = join(dirname(output_file), 'tmp.mustache.minimus.' + str(randint(1, 1e20)))
        for row in flank_rows:
            row['bam_path'] = bamfile
            row['outdir'] = tmp_outdir

        agents = threads
        with Pool(processes=agents) as pool:
            extensions = np.array(pool.map(extend, flank_rows))
        shell("rmdir %s" % tmp_outdir)

        flanks['extended'] = pd.DataFrame(flanks['consensus_seq'] != extensions)
        flanks['consensus_seq'] = extensions
        flanks['consensus_seq_length'] = pd.DataFrame(list(map(len, list(flanks['consensus_seq']))))

        flanks = flanks.loc[:, ['contig', 'pos', 'orient', 'softclip_count', 'runthrough_count',
                                'extended', 'consensus_seq_length', 'consensus_seq']]

        logger.info("Extended %d flanks using local assembly..." % sum(did_extend))

    if output_file:
        logger.info("Saving results to file %s" % output_file)
        flanks.to_csv(output_file, sep='\t', index=False)

    return flanks


@click.command()
@click.argument('flanksfile', type=click.Path(exists=True))
@click.argument('bamfile', type=click.Path(exists=True))
@click.option('--threads', '-t', default=1, help="The number of processors to run while finding flank extensions.")
@click.option('--output_file', '-o', default='mustache.extendflanks.tsv', help="The output file to save the results.")
def extendflanks(flanksfile, bamfile, threads, output_file=None):
    _extendflanks(flanksfile, bamfile, threads, output_file)


if __name__ == '__main__':
    extendflanks()
    #seqs = ['ACGCA', 'ACG', 'ACGC', 'ACGT', 'ACGTC', 'ACGTCA', 'ACGTCAT', 'ACGTCAG', 'ACGTCAT']
    #clusters = get_sequence_clusters(seqs)

    #mytrie = flanktrie.Trie()
    #for s in seqs:
    #    mytrie.add(s)

    #print(mytrie.traverse())
    #print(mytrie.calc_total_shared_words('ACGTCAG', 'ACGTCAT'))
    #print(mytrie.calc_total_unique_shared_words('ACGTCAG', 'ACGTCAT'))
    #print()

    #print(mytrie.traverse())
    #print("DELETING ACGTCAG")
    #mytrie.delete_word('ACG')
    #print(mytrie.traverse())