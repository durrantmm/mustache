import sys
import click
import pygogo as gogo
import pandas as pd
import numpy as np
from snakemake import shell
from mustache import fastatools, embosstools, fraggenescantools, hmmtools
from random import randint
from mustache.config import HMMDB


verbose=True
logger = gogo.Gogo(__name__, verbose=verbose).logger


# python mustache/findtransposase.py test/example_pairedflanks_tsv.txt 

# Run hmmsearch on the translated proteins using the ISEScan clusters.faa.hmm file
# Assign each sequence to an IS cluster if it has a hit
def _findtransposase(pairsfile, output_file=None, hmmdb=None):

	if not hmmdb:
		hmmdb = HMMDB
	
	flanks = pd.read_csv(pairsfile, sep='\t')
	
	logger.info("Writing flank pairs to fasta file...")
	tmp_flanks_fasta = '/tmp/mustache.findtransposase.' + str(randint(0, 1e100)) + '.fasta'
	fastatools.write_flanks_to_fasta(flanks, tmp_flanks_fasta)
	
	logger.info("Running FragGeneScan on flanks on %s..." % tmp_flanks_fasta)
	tmp_flanks_translate = '/tmp/mustache.findtransposase.' + str(randint(0, 1e100)) + '.translate'
	fraggenescantools.run_fraggenescan(tmp_flanks_fasta, tmp_flanks_translate)
	
	
	logger.info("Running hmmsearch on %s..." % tmp_flanks_translate)
	tmp_flanks_translate = tmp_flanks_translate + ".faa"
	tmp_hmm_results = '/tmp/mustache.findtransposase.' + str(randint(0, 1e100)) + '.hmmresults'
	hmmtools.run_hmmsearch(tmp_flanks_translate, tmp_hmm_results, hmmdb)
	
	hmm_results = hmmtools.process_hmm_results(tmp_hmm_results)
	
	final_results = flanks.join(hmm_results, how='left')
	
	#nohmmsearch = final_results.query("evalue_5p != evalue_5p & evalue_3p != evalue_3p")
		
	#logger.info("%d sequence pairs with HMM results...", final_results.shape[0] - nohmmsearch.shape[0])
	
	if output_file:
		logger.info("Saving results to file %s" % output_file)
		final_results.to_csv(output_file, sep='\t', index=False)
		
	shell('rm -f %s %s %s' % (tmp_flanks_fasta, tmp_flanks_translate, tmp_hmm_results))
	
	return final_results
    


@click.command()
@click.argument('pairsfile', type=click.Path(exists=True))
@click.option('--output_file', '-o', default='mustache.findtransposase.tsv', help="The output file to save the results.")
def findtransposase(pairsfile, output_file=None):
    _findtransposase(pairsfile, output_file)


if __name__ == '__main__':
    findtransposase()
