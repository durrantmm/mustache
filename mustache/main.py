import click
from mustache.alignbwa import alignbwa
from mustache.findflanks import findflanks, _findflanks
from mustache.pairflanks import pairflanks, _pairflanks
from mustache.blastpairs import blastpairs, _blastpairs
from mustache.blastpanisa import blastpanisa
from mustache.extendflanks import extendflanks, _extendflanks
from mustache.hmmsearchpairs import hmmsearchpairs, _hmmsearchpairs
from mustache.mergepairs import mergepairs
import pygogo as gogo
from os.path import isfile

verbose=True
logger = gogo.Gogo(__name__, verbose=verbose).logger

@click.group()
def cli():
    pass

@click.command()
@click.argument('bamfile', type=click.Path(exists=True))
@click.option('--output_prefix', '-o', default='mustache', help="The prefix to be used for all output files.")
@click.option('--min_softclip_length', '-minlen', default=4, help="For a softclipped site to be considered, there must be at least one softclipped read of this length.")
@click.option('--min_softclip_count', '-mincount', default=10, help="For a softclipped site to be considered, there must be at least this many softclipped reads at the site.")
@click.option('--min_alignment_quality', '-minq', default=20, help="For a read to be considered, it must meet this alignment quality cutoff.")
@click.option('--blastdb', '-blastdb')
@click.option('--checkexist/--no-checkexist', default=True)
def single_end_pipeline(bamfile, output_prefix, min_softclip_length, min_softclip_count, min_alignment_quality, blastdb,
                        checkexist):
    logger.info("BEGINNING FINDFLANKS...")
    findflanks_output_file = output_prefix + '.findflanks.tsv'
    if not (checkexist and isfile(findflanks_output_file)):
        checkexist = False
        _findflanks(bamfile, min_softclip_length, min_softclip_count, min_alignment_quality, findflanks_output_file)
    else:
        logger.info("OUTPUT FILE %s ALREADY EXISTS..." % findflanks_output_file)
    logger.info('')

    logger.info("BEGINNING PAIRFLANKS...")
    pairs_output_file = output_prefix + '.pairflanks.tsv'
    if not (checkexist and isfile(pairs_output_file)):
        checkexist = False
        _pairflanks(findflanks_output_file, pairs_output_file)
    else:
        logger.info("OUTPUT FILE %s ALREADY EXISTS..." % pairs_output_file)
    logger.info('')

    logger.info("BEGINNING BLASTPAIRS...")
    blastpairs_output_file = output_prefix + '.blastpairs.tsv'
    if not (checkexist and isfile(blastpairs_output_file)):
        checkexist = False
        _blastpairs(pairs_output_file, output_file=blastpairs_output_file)
    else:
        logger.info("OUTPUT FILE %s ALREADY EXISTS..." % blastpairs_output_file)
    logger.info('')

    logger.info("BEGINNING HMMSEARCHPAIRS...")
    hmmsearchpairs_output_file = output_prefix + '.hmmsearchpairs.tsv'
    if not (checkexist and isfile(hmmsearchpairs_output_file)):
        checkexist = False
        _hmmsearchpairs(blastpairs_output_file, output_file=hmmsearchpairs_output_file)
    else:
        logger.info("OUTPUT FILE %s ALREADY EXISTS..." % hmmsearchpairs_output_file)
    logger.info('')

    logger.info("MUSTACHE SINGLE-END PIPELINE COMPLETED SUCCESSFULLY.")


@click.command()
@click.argument('bamfile', type=click.Path(exists=True))
@click.option('--output_prefix', '-o', default='mustache', help="The prefix to be used for all output files.")
@click.option('--min_softclip_length', '-minlen', default=4, help="For a softclipped site to be considered, there must be at least one softclipped read of this length.")
@click.option('--min_softclip_count', '-mincount', default=10, help="For a softclipped site to be considered, there must be at least this many softclipped reads at the site.")
@click.option('--min_alignment_quality', '-minq', default=20, help="For a read to be considered, it must meet this alignment quality cutoff.")
@click.option('--blastdb', '-blastdb')
@click.option('--checkexist/--no-checkexist', default=True)
def paired_end_pipeline(bamfile, output_prefix, min_softclip_length, min_softclip_count, min_alignment_quality, blastdb,
                        checkexist):

    logger.info("BEGINNING FINDFLANKS...")
    findflanks_output_file = output_prefix + '.findflanks.tsv'
    if not (checkexist and isfile(findflanks_output_file)):
        checkexist = False
        _findflanks(bamfile, min_softclip_length, min_softclip_count, min_alignment_quality, findflanks_output_file)
    else:
        logger.info("OUTPUT FILE %s ALREADY EXISTS..." % findflanks_output_file)
    logger.info('')

    logger.info("BEGINNING EXTENDFLANKS...")
    extendflanks_output_file = output_prefix + '.extendflanks.tsv'
    if not (checkexist and isfile(extendflanks_output_file)):
        checkexist = False
        _extendflanks(findflanks_output_file, bamfile, extendflanks_output_file)
    else:
        logger.info("OUTPUT FILE %s ALREADY EXISTS..." % extendflanks_output_file)
    logger.info('')

    logger.info("BEGINNING PAIRFLANKS...")
    pairs_output_file = output_prefix + '.pairflanks.tsv'
    if not (checkexist and isfile(pairs_output_file)):
        checkexist = False
        _pairflanks(extendflanks_output_file, pairs_output_file)
    else:
        logger.info("OUTPUT FILE %s ALREADY EXISTS..." % pairs_output_file)
    logger.info('')

    logger.info("BEGINNING BLASTPAIRS...")
    blastpairs_output_file = output_prefix + '.blastpairs.tsv'
    if not (checkexist and isfile(blastpairs_output_file)):
        checkexist = False
        _blastpairs(pairs_output_file, output_file=blastpairs_output_file)
    else:
        logger.info("OUTPUT FILE %s ALREADY EXISTS..." % blastpairs_output_file)
    logger.info('')

    logger.info("BEGINNING HMMSEARCHPAIRS...")
    hmmsearchpairs_output_file = output_prefix + '.hmmsearchpairs.tsv'
    if not (checkexist and isfile(hmmsearchpairs_output_file)):
        checkexist = False
        _hmmsearchpairs(blastpairs_output_file, output_file=hmmsearchpairs_output_file )
    else:
        logger.info("OUTPUT FILE %s ALREADY EXISTS..." % hmmsearchpairs_output_file )
    logger.info('')


    logger.info("MUSTACHE PAIRED-END PIPELINE COMPLETED SUCCESSFULLY.")


cli.add_command(alignbwa)
cli.add_command(findflanks)
cli.add_command(pairflanks)
cli.add_command(extendflanks)
cli.add_command(blastpairs)
cli.add_command(hmmsearchpairs)
cli.add_command(blastpanisa)
cli.add_command(single_end_pipeline)
cli.add_command(mergepairs)
cli.add_command(paired_end_pipeline)


if __name__ == '__main__':

    cli()