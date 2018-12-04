import warnings
warnings.filterwarnings("ignore")
import click
from mustache.formatbam import formatbam
from mustache.findflanks import findflanks, _findflanks
from mustache.pairflanks import pairflanks, _pairflanks
from mustache.extendpairs import extendpairs, _extendpairs
from mustache.recall import recall
from mustache.inferseqassembly import inferseq_assembly
from mustache.inferseqreference import inferseq_reference
from mustache.inferseqdatabase import inferseq_database
from mustache.inferseqoverlap import inferseq_overlap

import pygogo as gogo
from os.path import isfile
from os.path import basename, dirname
from os import makedirs

verbose=True
logger = gogo.Gogo(__name__, verbose=verbose).logger


@click.group()
def cli():
    pass


@click.command()
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


cli.add_command(formatbam)
cli.add_command(findflanks)
cli.add_command(pairflanks)
cli.add_command(pipeline)
cli.add_command(recall)
cli.add_command(inferseq_assembly)
cli.add_command(inferseq_reference)
cli.add_command(inferseq_database)
cli.add_command(inferseq_overlap)
cli.add_command(extendpairs)


if __name__ == '__main__':

    cli()