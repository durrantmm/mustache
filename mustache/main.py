import warnings
warnings.filterwarnings("ignore")
import click

from mustache.findflanks import _findflanks
from mustache.pairflanks import _pairflanks
from mustache.extendpairs import _extendpairs
from mustache.inferseqassembly import _inferseq_assembly
from mustache.inferseqoverlap import _inferseq_overlap
from mustache.inferseqreference import _inferseq_reference
from mustache.inferseqdatabase import _inferseq_database
from mustache.formatbam import _formatbam
from mustache.recall import _recall
from mustache.help import CustomHelp

import pygogo as gogo
from os.path import isfile
from os.path import basename, dirname
from os import makedirs

verbose=True
logger = gogo.Gogo(__name__, verbose=verbose).logger


@click.group(cls=CustomHelp)
def cli():
    """Command-line tools to identify mobile element insertions from short-read sequencing data."""
    pass

@cli.command(short_help='Find insertion sites and reconstruct flanks of inserted sequence', help_priority=1)
@click.argument('bamfile', type=click.Path(exists=True))
@click.option('--output_file', '-o', default='mustache.findflanks.tsv', help="The output file to save the results.")
@click.option('--min_softclip_length', '-minlen', default=4, help="For a softclipped site to be considered, there must be at least one softclipped read of this length.")
@click.option('--min_softclip_count', '-mincount', default=4, help="For a softclipped site to be considered, there must be at least this many softclipped reads at the site.")
@click.option('--min_count_consensus', '-mcc', default=2, help="When building the consensus sequence, stop building consensus if read count drops below this cutoff.")
@click.option('--min_alignment_quality', '-minq', default=20, help="For a read to be considered, it must meet this alignment quality cutoff.")
@click.option('--min_softclip_ratio', '-minratio', default=0.1, help="For a softclipped site to be considered, the proportion of softclipped sites must not fall below this value. Below 0 and 0.1.")
@click.option('--min_alignment_inner_length', '-minial', default=21, help="If a read is softclipped on both ends, the aligned portion must be at least this long. Ideally, set this equal to 1 + maximum direct repeat length.")
def findflanks(bamfile, min_softclip_length, min_softclip_count, min_count_consensus, min_alignment_quality, min_softclip_ratio, min_alignment_inner_length, output_file=None):
    """A click access point for the findflanks module. This is used for creating the command line interface."""

    _findflanks(bamfile, min_softclip_length, min_softclip_count, min_count_consensus, min_alignment_quality, min_softclip_ratio,
                min_alignment_inner_length, output_file)


@cli.command(short_help="Pair identified flanks with each other to represent 5' and 3' ends of inserted sequence.", help_priority=2)
@click.argument('flanksfile', type=click.Path(exists=True))
@click.argument('bamfile', type=click.Path(exists=True))
@click.argument('genome', type=click.Path(exists=True))
@click.option('--max_direct_repeat_length', '-maxdr', default=20, help="The maximum length of a direct repeat to consider a pair.")
@click.option('--output_file', '-o', default='mustache.pairflanks.tsv', help="The output file to save the results.")
def pairflanks(flanksfile, bamfile, genome, max_direct_repeat_length, output_file=None):
    _pairflanks(flanksfile, bamfile, genome, max_direct_repeat_length, output_file)


@cli.command(help_priority=3)
@click.argument('pairsfile', type=click.Path(exists=True))
@click.argument('bamfile', type=click.Path(exists=True))
@click.option('--threads', '-t', default=1, help="The number of processors to run while finding flank extensions.")
@click.option('--output_file', '-o', default='mustache.extendpairs.tsv', help="The output file to save the results.")
def extendpairs(pairsfile, bamfile, threads, output_file=None):
    """
    Experimental. Extends the consensus flanks using a local assembly of paired end reads.
    BAM file must be processed using the 'formatbam' command first.
    Requires an installation of the AMOS sequence assembly software: http://amos.sourceforge.net/wiki/index.php/AMOS
    """
    _extendpairs(pairsfile, bamfile, threads, output_file)


@cli.command(short_help='Infers the identity of an inserted sequence by aligning flank pairs to an assembled genome.', help_priority=4)
@click.argument('pairsfile', type=click.Path(exists=True))
@click.argument('bamfile', type=click.Path(exists=True))
@click.argument('inferseq_assembly', type=click.Path(exists=True))
@click.argument('inferseq_reference', type=click.Path(exists=True))
@click.option('--min_perc_identity', '-minident', default=0.95, help="Only consider matches with a percentage identity above this threshold")
@click.option('--max_internal_softclip_prop', '-maxclip', default=0.05, help="Do not consider matches with internal softclipped ends exceeding this proportion of the total read")
@click.option('--max_inferseq_size', '-maxsize', default=500000, help="Do not consider inferred sequences over this size.")
@click.option('--min_inferseq_size', '-minsize', default=1, help="Do not consider inferred sequences below this size.")
@click.option('--keep-intermediate/--no-keep-intermediate', default=False, help="Keep intermediate files.")
@click.option('--output_file', '-o', default='mustache.inferseq_assembly.tsv', help="The output file to save the results.")
def inferseq_assembly(pairsfile, bamfile, inferseq_assembly, inferseq_reference, min_perc_identity,
                      max_internal_softclip_prop, max_inferseq_size, min_inferseq_size, keep_intermediate, output_file=None):
    """Infers the identity of an inserted sequence by aligning flank pairs to an assembled genome."""

    _inferseq_assembly(pairsfile, bamfile, inferseq_assembly, inferseq_reference, min_perc_identity,
                       max_internal_softclip_prop, max_inferseq_size, min_inferseq_size, keep_intermediate, output_file)


@cli.command(short_help='Infers the identity of an inserted sequence by aligning flank pairs to a reference genome. Ideal for re-sequencing experiments where evolved strains are closely related to the reference genome used.',
             help_priority = 5)
@click.argument('pairsfile', type=click.Path(exists=True))
@click.argument('inferseq_reference', type=click.Path(exists=True))
@click.option('--min_perc_identity', '-minident', default=0.95, help="Only consider matches with a percentage identity above this threshold")
@click.option('--max_internal_softclip_prop', '-maxclip', default=0.05, help="Do not consider matches with internal softclipped ends exceeding this proportion of the total read")
@click.option('--max_inferseq_size', '-maxsize', default=500000, help="Do not consider inferred sequences over this size.")
@click.option('--min_inferseq_size', '-minsize', default=1, help="Do not consider inferred sequences below this size.")
@click.option('--output_file', '-o', default='mustache.inferseq_database.tsv', help="The output file to save the results.")
@click.option('--keep-intermediate/--no-keep-intermediate', default=False, help="Keep intermediate files.")
@click.option('--output_file', '-o', default='mustache.inferseq_reference.tsv', help="The output file to save the results.")
def inferseq_reference(pairsfile, inferseq_reference, min_perc_identity, max_internal_softclip_prop,
                       max_inferseq_size, min_inferseq_size, keep_intermediate, output_file=None):
    """
    Infers the identity of an inserted sequence by aligning flank pairs to a reference genome.
    Ideal for re-sequencing experiments where evolved strains are closely related to the reference genome used.
    """

    _inferseq_reference(pairsfile, inferseq_reference, min_perc_identity, max_internal_softclip_prop,
                        max_inferseq_size, min_inferseq_size, keep_intermediate, output_file)


@cli.command(short_help='Infers the identity of an inserted sequence by checking if they overlap with one another. Only identifies an inserted sequence if the consensus flanks are long enough to span the entire insertion.',
             help_priority=6)
@click.argument('pairsfile', type=click.Path(exists=True))
@click.option('--min_overlap_score', '-minscore', default=10, help="The minimum overlap score to keep inferred sequence.")
@click.option('--min_overlap_perc_identity', '-min_overlap_perc_identity', default=0.9, help="The minimum overlap percent identity to keep inferred sequence")
@click.option('--output_file', '-o', default='mustache.inferseq_overlap.tsv', help="The output file to save the results.")
def inferseq_overlap(pairsfile, min_overlap_score, min_overlap_perc_identity, output_file=None):
    """
    Infers the identity of an inserted sequence by checking if they overlap with one another.
    Only identifies an inserted sequence if the consensus flanks are long enough to span the entire insertion.
    """

    _inferseq_overlap(pairsfile, min_overlap_score, min_overlap_perc_identity, output_file)


@cli.command(short_help='Infers the identity of an inserted sequence by aligning flank pairs to an database of known inserted elements.', help_priority=7)
@click.argument('pairsfile', type=click.Path(exists=True))
@click.argument('fasta_database', type=click.Path(exists=True))
@click.option('--min_perc_identity', '-minident', default=0.95, help="Only consider matches with a percentage identity above this threshold")
@click.option('--max_internal_softclip_prop', '-maxclip', default=0.05, help="Do not consider matches with internal softclipped ends exceeding this proportion of the total read")
@click.option('--max_edge_distance', '-maxedgedist', default=10, help="Reads must align within this number of bases from the edge of an element to be considered.")
@click.option('--output_file', '-o', default='mustache.inferseq_database.tsv', help="The output file to save the results.")
@click.option('--keep-intermediate/--no-keep-intermediate', default=False, help="Keep intermediate files.")
def inferseq_database(pairsfile, fasta_database, min_perc_identity,  max_internal_softclip_prop, max_edge_distance, output_file=None, keep_intermediate=False):
    """Infers the identity of an inserted sequence by aligning flank pairs to an database of known inserted elements."""

    _inferseq_database(pairsfile, fasta_database, min_perc_identity, max_internal_softclip_prop, max_edge_distance, output_file, keep_intermediate)


@cli.command(short_help="Formats a BAM file for use with mustache. Usually not necessary, unless using the experiment extendpairs command.", help_priority=8)
@click.argument('in_sam', type=click.Path(exists=True))
@click.argument('out_bam')
@click.option('--single-end', is_flag=True, default=False, help="Add this flag for single-end files")
@click.option('--keep-tmp-files', is_flag=True, default=False, help="Add this flag if you want to keep intermediate temporary files")
def formatbam(in_sam, out_bam, single_end, keep_tmp_files):

    _formatbam(in_sam, out_bam, single_end, keep_tmp_files)



@cli.command(short_help='Recall softclip counts and runthrough counts from BAM file at specified pairflank insertions.', help_priority=9)
@click.argument('pairsfile', type=click.Path(exists=True))
@click.argument('bamfile', type=click.Path(exists=True))
@click.option('--min_alignment_quality', '-minq', default=20, help="For a read to be considered, it must meet this alignment quality cutoff.")
@click.option('--min_alignment_inner_length', '-minial', default=21, help="If a read is softclipped on both ends, the aligned portion must be at least this long. Ideally, set this equal to 1 + maximum direct repeat length.")
@click.option('--output_file', '-o', default='mustache.recall.tsv', help="The output file to save results to...")
def recall(pairsfile, bamfile, min_alignment_quality, min_alignment_inner_length, output_file):
    _recall(pairsfile, bamfile, min_alignment_quality, min_alignment_inner_length, output_file)



@cli.command(help_priority=10)
@click.argument('bamfile', type=click.Path(exists=True))
@click.argument('genome', type=click.Path(exists=True))
@click.option('--output_prefix', '-o', default=None, help="The prefix to be used for all output files.")
@click.option('--threads', '-t', default=1, help="Number of threads to use to extend flanks in paired sequences.")
@click.option('--min_softclip_length', '-minlen', default=4, help="For a softclipped site to be considered, there must be at least one softclipped read of this length.")
@click.option('--min_softclip_count', '-mincount', default=4, help="For a softclipped site to be considered, there must be at least this many softclipped reads at the site.")
@click.option('--min_alignment_quality', '-minq', default=20, help="For a read to be considered, it must meet this alignment quality cutoff.")
@click.option('--min_softclip_ratio', '-minratio', default=0.1, help="For a softclipped site to be considered, the ratio of softclipped sites runthrough sites must not fall below this value.")
@click.option('--max_direct_repeat_length', '-maxdr', default=20, help="The maximum length of a direct repeat to consider a pair.")
@click.option('--blastdb', '-blastdb')
@click.option('--checkexist/--no-checkexist', default=True)
@click.option('--annotate/--no-annotate', default=False)
def pipeline(bamfile, genome, output_prefix, threads, min_softclip_length, min_softclip_count,
             min_alignment_quality, min_softclip_ratio, max_direct_repeat_length, blastdb, checkexist, annotate):

    min_alignment_inner_length = max_direct_repeat_length + 1

    if output_prefix:
        makedirs(dirname(output_prefix), exist_ok=True)

    if not output_prefix:
        output_prefix = '.'.join(['mustache', basename(bamfile).split('.')[0]])


    logger.info("BEGINNING FINDFLANKS...")
    findflanks_output_file = output_prefix + '.findflanks.tsv'
    if not (checkexist and isfile(findflanks_output_file)):
        checkexist = False
        _findflanks(bamfile, min_softclip_length, min_softclip_count, min_alignment_quality,
                    min_softclip_ratio, min_alignment_inner_length, findflanks_output_file)
    else:
        logger.info("OUTPUT FILE %s ALREADY EXISTS..." % findflanks_output_file)
    logger.info('')

    logger.info("BEGINNING PAIRFLANKS...")
    pairs_output_file = output_prefix + '.pairflanks.tsv'
    if not (checkexist and isfile(pairs_output_file)):
        checkexist = False
        _pairflanks(findflanks_output_file, bamfile, genome, max_direct_repeat_length, pairs_output_file)
    else:
        logger.info("OUTPUT FILE %s ALREADY EXISTS..." % pairs_output_file)
    logger.info('')

    logger.info("BEGINNING EXTENDPAIRS...")
    extendpairs_output_file = output_prefix + '.extendpairs.tsv'
    if not (checkexist and isfile(extendpairs_output_file)):
        checkexist = False
        _extendpairs(pairs_output_file, bamfile, threads, extendpairs_output_file)
    else:
        logger.info("OUTPUT FILE %s ALREADY EXISTS..." % extendpairs_output_file)
    logger.info('')

    logger.info("MUSTACHE PIPELINE COMPLETED SUCCESSFULLY.")



cli.commands['findflanks'].help_priority = 1

if __name__ == '__main__':

    cli()