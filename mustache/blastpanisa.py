import warnings
warnings.filterwarnings("ignore")
import sys
import click
import pygogo as gogo
import pandas as pd
import numpy as np
from snakemake import shell
from mustache import fastatools, embosstools, blasttools
from random import randint
from mustache.config import BLASTDB
from os.path import join

verbose=True
logger = gogo.Gogo(__name__, verbose=verbose).logger


def _blastpanisa(panisafile, blastdb=None, output_file=None):

    if not blastdb:
        blastdb = BLASTDB

    panisa = pd.read_csv(panisafile, sep='\t')

    logger.info("Writing flank pairs to fasta file...")
    tmp_flanks_fasta = join(TMPDIR, 'mustache.blastpanisa.' + str(randint(0, 1e20)) + '.fasta')
    fastatools.write_panisa_flanks_to_fasta(panisa, tmp_flanks_fasta)

    logger.info("Blasting panISa flank pairs against database found at %s" % blastdb)
    tmp_flanks_blast = join(TMPDIR, 'mustache.blastpanisa.' + str(randint(0, 1e20)) + '.blast.txt')
    blasttools.blast_fasta(tmp_flanks_fasta, blastdb, tmp_flanks_blast)

    logger.info("Processing blast results %s" % blastdb)
    blast_results = blasttools.process_blast_results(tmp_flanks_blast)

    final_results = panisa.join(blast_results, how='left')

    noblast = final_results.query("evalue_5p != evalue_5p & evalue_3p != evalue_3p")

    logger.info("%d sequence pairs blasted to an IS...", final_results.shape[0] - noblast.shape[0])
    if output_file:
        logger.info("Saving results to file %s" % output_file)
        final_results.to_csv(output_file, sep='\t', index=False)

    shell('rm -f %s %s' % (tmp_flanks_fasta, tmp_flanks_blast))

    return final_results




@click.command()
@click.argument('panisafile', type=click.Path(exists=True))
@click.option('--blastdb', '-db')
@click.option('--output_file', '-o', default='mustache.blastpairs.tsv', help="The output file to save the results.")
def blastpanisa(panisafile, blastdb=None, output_file=None):
    _blastpanisa(panisafile, blastdb, output_file)


if __name__ == '__main__':
    blastpanisa()