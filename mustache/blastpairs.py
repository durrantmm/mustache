import warnings
warnings.filterwarnings("ignore")
import sys
import click
import pygogo as gogo
import pandas as pd
import numpy as np
from snakemake import shell
from mustache import fastatools, blasttools
from random import randint
from mustache.config import BLASTDB, TMPDIR
from os.path import basename, join

verbose=True
logger = gogo.Gogo(__name__, verbose=verbose).logger


def _blastpairs(pairsfile, blastdb=None, output_file=None):

    if not blastdb:
        blastdb = BLASTDB

    pairs = pd.read_csv(pairsfile, sep='\t')

    if pairs.shape[0] == 0:
        logger.info("No pairs found in the input file...")

        if output_file:
            logger.info("Saving results to file %s" % output_file)
            pairs.to_csv(output_file, sep='\t', index=False)

        return pairs

    else:
        logger.info("Writing flank pairs to fasta file...")
        tmp_pairs_fasta = join(TMPDIR, 'mustache.blastpairs.' + str(randint(0, 1e20)) + '.fasta')
        fastatools.write_flanks_to_fasta(pairs, tmp_pairs_fasta)

        logger.info("Blasting flank pairs against database found at %s" % blastdb)
        tmp_pairs_blast = join(TMPDIR, 'mustache.blastpairs.' + str(randint(0, 1e20)) + '.blast.txt')
        blasttools.blast_fasta(tmp_pairs_fasta, blastdb, tmp_pairs_blast)

        logger.info("Processing blast results %s" % blastdb)
        blast_results = blasttools.process_blast_results(tmp_pairs_blast)

        final_results = pairs.join(blast_results, how='left')

        noblast = final_results.query("blast_evalue_5p != blast_evalue_5p & blast_evalue_3p != blast_evalue_3p")

        logger.info("%d sequence pairs blasted to an IS...", final_results.shape[0] - noblast.shape[0])
        if output_file:
            logger.info("Saving results to file %s" % output_file)
            final_results.to_csv(output_file, sep='\t', index=False)

        shell('rm -f %s %s' % (tmp_pairs_fasta, tmp_pairs_blast))

    return final_results




@click.command()
@click.argument('pairsfile', type=click.Path(exists=True))
@click.option('--blastdb', '-db')
@click.option('--output_file', '-o', default=None, help="The output file to save the results.")
def blastpairs(pairsfile, blastdb=None, output_file=None):
    if not output_file:
        output_file = '.'.join(['mustache', basename(pairsfile).split('.')[0], 'blastpairs.tsv'])
    _blastpairs(pairsfile, blastdb, output_file)


if __name__ == '__main__':
    blastpairs()