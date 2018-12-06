import warnings
warnings.filterwarnings("ignore")
import sys
import click
from mustache import pysamtools
import pandas as pd
import pysam
import pygogo as gogo

verbose=True
logger = gogo.Gogo(__name__, verbose=verbose).logger

def _recall(pairsfile, bamfile, min_alignment_quality, min_alignment_inner_length, output_file):
    pairs = pd.read_csv(pairsfile, sep='\t')
    bam = pysam.AlignmentFile(bamfile)

    logger.info("Getting runthrough counts at all unique sites...")
    pairs = count_runthrough_reads(pairs, bam, min_alignment_quality, min_alignment_inner_length)
    logger.info("Getting softclipped counts at all unique sites...")
    pairs = count_softclipped_reads(pairs, bam, min_alignment_quality, min_alignment_inner_length)

    pairs.to_csv(output_file, sep='\t', index=False)


def count_runthrough_reads(pairs, bam, min_qual, min_alignment_inner_length):

    threeprime_sites = pairs[['contig', 'pos_3p']].drop_duplicates()
    threeprime_sites.columns = ['contig', 'pos']
    fiveprime_sites = pairs[['contig', 'pos_5p']].drop_duplicates()
    fiveprime_sites.columns = ['contig', 'pos']
    unique_sites = pd.concat([threeprime_sites, fiveprime_sites]).drop_duplicates()

    runthrough_counts = []
    for index, row in unique_sites.iterrows():
        runthrough_counts.append(pysamtools.count_runthrough_reads(bam, row['contig'], row['pos'], min_qual,
                                                                   min_alignment_inner_length))

    unique_sites.loc[:, 'runthrough_count'] = runthrough_counts

    threeprime_unique_sites = unique_sites
    threeprime_unique_sites.columns = ['contig', 'pos_3p', 'runthrough_count_3p']
    outpairs = pairs.merge(threeprime_unique_sites, how='left')

    fiveprime_unique_sites = unique_sites
    fiveprime_unique_sites.columns = ['contig', 'pos_5p', 'runthrough_count_5p']
    outpairs = outpairs.merge(fiveprime_unique_sites, how='left')

    return outpairs

def count_softclipped_reads(pairs, bam, min_qual, min_alignment_inner_length):

    threeprime_sites = pairs[['contig', 'pos_3p']].drop_duplicates()
    threeprime_sites.columns = ['contig', 'pos']
    fiveprime_sites = pairs[['contig', 'pos_5p']].drop_duplicates()
    fiveprime_sites.columns = ['contig', 'pos']
    unique_sites = pd.concat([threeprime_sites, fiveprime_sites]).drop_duplicates()

    softclip_counts = []
    for index, row in unique_sites.iterrows():
        softclip_counts.append(pysamtools.count_softclipped_reads(bam, row['contig'], row['pos'], min_qual,
                                                                   min_alignment_inner_length))

    unique_sites.loc[:, 'softclip_count'] = softclip_counts

    threeprime_unique_sites = unique_sites
    threeprime_unique_sites.columns = ['contig', 'pos_3p', 'softclip_count_3p']
    outpairs = pairs.merge(threeprime_unique_sites, how='left')

    fiveprime_unique_sites = unique_sites
    fiveprime_unique_sites.columns = ['contig', 'pos_5p', 'softclip_count_5p']
    outpairs = outpairs.merge(fiveprime_unique_sites, how='left')

    return outpairs

if __name__ == '__main__':
    recall()