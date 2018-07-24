import sys
import click
import pygogo as gogo
import pandas as pd
import numpy as np
from snakemake import shell
from mustache import fastatools, embosstools, fraggenescantools, hmmtools
from random import randint
from mustache.config import HMMDB
from os.path import isfile

verbose = True
logger = gogo.Gogo(__name__, verbose=verbose).logger


# python mustache/hmmsearchpairs.py test/example_pairedflanks_tsv.txt 

# Run hmmsearch on the translated proteins using the ISEScan clusters.faa.hmm file
# Assign each sequence to an IS cluster if it has a hit
def _hmmsearchpairs(pairsfile, output_file=None, hmmdb=None):
    if not hmmdb:
        hmmdb = HMMDB

    pairs = pd.read_csv(pairsfile, sep='\t')

    if pairs.shape[0] == 0:
        logger.info("No pairs found in the input file...")

        if output_file:
            logger.info("Saving results to file %s" % output_file)
            pairs.to_csv(output_file, sep='\t', index=False)

        return pairs

    logger.info("Writing flank pairs to fasta file...")
    tmp_pairs_fasta = '/tmp/mustache.hmmsearchpairs.' + str(randint(0, 1e100)) + '.fasta'
    fastatools.write_flanks_to_fasta(pairs, tmp_pairs_fasta)

    logger.info("Running FragGeneScan on pairs on %s..." % tmp_pairs_fasta)
    tmp_pairs_translate = '/tmp/mustache.hmmsearchpairs.' + str(randint(0, 1e100)) + '.translate'
    fraggenescantools.run_fraggenescan(tmp_pairs_fasta, tmp_pairs_translate)

    has_fraggene_output = True
    if fraggenescantools.count_fraggenescan_seqs(tmp_pairs_translate+'.faa') == 0:
        has_fraggene_output = False

    logger.info("Running hmmsearch on %s..." % tmp_pairs_translate)
    tmp_pairs_translate = tmp_pairs_translate + ".faa"
    tmp_hmm_results = '/tmp/mustache.hmmsearchpairs.' + str(randint(0, 1e100)) + '.hmmresults'

    if has_fraggene_output:
        hmmtools.run_hmmsearch(tmp_pairs_translate, tmp_hmm_results, hmmdb)

    hmm_results = hmmtools.process_hmm_results(tmp_hmm_results)

    final_results = pairs.join(hmm_results, how='left')

    # nohmmsearch = final_results.query("evalue_5p != evalue_5p & evalue_3p != evalue_3p")

    # logger.info("%d sequence pairs with HMM results...", final_results.shape[0] - nohmmsearch.shape[0])

    if output_file:
        logger.info("Saving results to file %s" % output_file)
        final_results.to_csv(output_file, sep='\t', index=False)

    shell('rm -f %s %s %s' % (tmp_pairs_fasta, tmp_pairs_translate, tmp_hmm_results))

    return final_results


@click.command()
@click.argument('pairsfile', type=click.Path(exists=True))
@click.option('--output_file', '-o', default='mustache.hmmsearchpairs.tsv', help="The output file to save the results.")
def hmmsearchpairs(pairsfile, output_file=None):
    _hmmsearchpairs(pairsfile, output_file)


if __name__ == '__main__':
    hmmsearchpairs()
