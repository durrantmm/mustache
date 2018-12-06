import warnings
warnings.filterwarnings("ignore")
import sys
import click
from mustache import misc, sctools, flanktrie, pysamtools
import pysam
from collections import defaultdict
import pandas as pd
import pygogo as gogo
from jellyfish import levenshtein_distance
import numpy as np
from scipy.sparse.csgraph import connected_components
import itertools
from os.path import basename
verbose=True
logger = gogo.Gogo(__name__, verbose=verbose).logger
from scipy.stats import binom_test

import warnings
warnings.filterwarnings("ignore", message="numpy.dtype size changed")
warnings.filterwarnings("ignore", message="numpy.ufunc size changed")

def get_softclipped_sites(bam_file, min_softclip_length=4, min_softclip_count=10, min_alignment_quality=20,
                          min_alignment_inner_length=21, min_softclip_ratio=0.1):

    column_names = ["contig", "pos", "orient", "softclip_count", "meets_minlength"]
    column_index_key = dict(list(zip(column_names, range(len(column_names)))))

    soft_clips = dict()

    read_count = 0
    for read in bam_file.fetch():
        read_count += 1
        if read_count % 100000 == 0:
            logger.info("\t\tAfter checking %d reads, %d softclipped sites found..." % (read_count, len(soft_clips)))

        if read.mapping_quality < min_alignment_quality:
            continue

        if not sctools.read_meets_min_alignment_inner_length(read, min_alignment_inner_length):
            continue

        if sctools.is_right_softclipped_strict(read):

            meets_minlength = False
            if sctools.right_softclip_length(read) >= min_softclip_length:
                meets_minlength = True

            contig, pos = sctools.right_softclipped_site_lenient(read)
            orient = 'R'
            index = (contig, pos, orient)

            if index in soft_clips:
                soft_clips[index][column_index_key['softclip_count']] += 1
                if soft_clips[index][column_index_key['meets_minlength']] == False and meets_minlength:
                    soft_clips[index][column_index_key['meets_minlength']] = meets_minlength
            else:
                new_site = [contig, pos, orient, 1, meets_minlength]
                soft_clips[index] = new_site

        if sctools.is_left_softclipped_strict(read):

            meets_minlength = False
            if sctools.left_softclip_length(read) >= min_softclip_length:
                meets_minlength = True

            contig, pos = sctools.left_softclipped_site_lenient(read)
            orient = 'L'
            index = (contig, pos, orient)

            if index in soft_clips:
                soft_clips[index][column_index_key['softclip_count']] += 1
                if soft_clips[index][column_index_key['meets_minlength']] == False and meets_minlength:
                    soft_clips[index][column_index_key['meets_minlength']] = meets_minlength
            else:
                new_site = [contig, pos, orient, 1, meets_minlength]
                soft_clips[index] = new_site

    logger.info("Applying initial filters to softclipped sites...")
    soft_clips = pd.DataFrame.from_dict(soft_clips, orient='index', columns=column_names).sort_index()

    soft_clips.query("meets_minlength == True & softclip_count > @min_softclip_count", inplace=True)

    soft_clips = filter_lone_flanks(soft_clips)
    logger.info("\tAfter removing lone softclips, filtering by minimum softclip length of %d, "
                "and filtering by a minimum softclip count of %d, "
                "%d sites remaining..." % (min_softclip_length, min_softclip_count, soft_clips.shape[0]))


    if soft_clips.shape[0] == 0:
        return soft_clips

    logger.info("Getting runthrough reads at remaining sites...")
    soft_clips_runthrough = count_runthrough_reads(soft_clips, bam_file, min_alignment_quality, min_alignment_inner_length)

    soft_clips_runthrough.index = soft_clips.index
    soft_clips = soft_clips_runthrough

    soft_clips.loc[:,'ratio'] = soft_clips.apply(
        lambda x: x['softclip_count'] / (x['softclip_count']+x['runthrough_count']), axis=1)

    soft_clips.query("ratio >= @min_softclip_ratio", inplace=True)
    soft_clips = soft_clips.loc[:, ['contig', 'pos', 'orient', 'softclip_count', 'runthrough_count', 'meets_minlength']]


    if soft_clips.shape[0] == 0:
        return soft_clips

    logger.info("\tAfter filtering by a count ratio cutoff 0f %f, %d sites remaining..." % (min_softclip_ratio, soft_clips.shape[0]))
    return soft_clips


def determine_flank_sequence(softclipped_sites, bam_file, min_softclip_length, min_softclip_count, min_count_consensus,
                             min_alignment_quality, min_alignment_inner_length, min_softclip_ratio):

    column_names = ["contig", "pos", "orient", "softclip_count", "consensus_seq_length", "consensus_seq"]

    runthrough_counts = softclipped_sites.loc[:,['contig', 'pos', 'runthrough_count']].drop_duplicates().reset_index(drop=True)

    final_flank_seqs = dict()
    site_count = 0
    for index, row in softclipped_sites.iterrows():
        site_count += 1
        if site_count % 1000 == 0:
            logger.info("\tAfter processing %d sites, %d unique softclipped flanks found..." % (site_count, len(final_flank_seqs)))

        contig, pos, orient = index


        softclipped_seqs, softclipped_qualities = pysamtools.get_softclipped_seqs_qualities(
            bam_file, contig, pos, orient, min_alignment_quality, min_alignment_inner_length)

        flanktrie = get_flank_sequence_trie(softclipped_seqs, softclipped_qualities)
        flanktrie_traversal = flanktrie.traverse_all()
        seq_dict = {seq[0]:[seq[1], seq[2]] for seq in flanktrie_traversal}

        seq_clusters = get_sequence_clusters(seq_dict)

        consensus_seqs = get_cluster_consensus_seqs(seq_clusters, min_count_consensus)

        if orient == 'L':
            consensus_seqs = [seq[::-1] for seq in consensus_seqs]

        merged_counts = merge_cluster_counts(seq_clusters, flanktrie)


        for seq, count in zip(consensus_seqs, merged_counts):
            final_flank_seqs[len(final_flank_seqs)+1] = [contig, pos, orient, count, len(seq), seq]

    final_flank_seqs = pd.DataFrame.from_dict(final_flank_seqs, orient='index', columns=column_names).sort_index()
    final_flank_seqs.query('consensus_seq_length >= @min_softclip_length & softclip_count >= @min_softclip_count', inplace=True)

    final_flank_seqs = pd.merge(final_flank_seqs, runthrough_counts, how='left')

    final_flank_seqs.loc[:, 'ratio'] = final_flank_seqs.apply(
        lambda x: x['softclip_count'] / (x['softclip_count'] + x['runthrough_count']), axis=1)

    final_flank_seqs.query("ratio >= @min_softclip_ratio", inplace=True)

    logger.info("\tAfter filtering by minimum softclip length of %d, "
                "a minimum softclip count of %d, and a count ratio cutoff 0f %f, "
                "%d sites remaining..." % (
                min_softclip_length, min_softclip_count, min_softclip_ratio, final_flank_seqs.shape[0]))

    return final_flank_seqs


def get_flank_sequence_trie(softclipped_seqs, softclipped_base_qual):

    mytrie = flanktrie.Trie()

    for i in range(len(softclipped_seqs)):
        seq, qual = softclipped_seqs[i], softclipped_base_qual[i]
        mytrie.add(seq.upper(), qual)

    return mytrie


def get_sequence_clusters(sequence_dict, perc_similarity=0.75):

    if len(sequence_dict) == 1:
        return [sequence_dict]

    sequences = list(sequence_dict.keys())

    no_match = set(sequences)
    has_match = set()
    pairwise_matches = set()
    for i, j in itertools.combinations(range(len(sequences)), r=2):

        seq1, seq2 = sequences[i], sequences [j]
        minlength = min([len(seq1), len(seq2)])
        dist = levenshtein_distance(seq1[:minlength], seq2[:minlength]     )
        similarity = 1 - (dist / minlength)

        if similarity >= perc_similarity:
            has_match.add(i)
            has_match.add(j)
            if seq1 in no_match: no_match.remove(seq1)
            if seq2 in no_match: no_match.remove(seq2)
            pairwise_matches.add((i, j))

    pair_order = sorted(list(has_match))
    adjacency_matrix = pd.DataFrame(np.zeros((len(pair_order), len(pair_order))), index=pair_order, columns=pair_order)

    for pair in pairwise_matches:
        adjacency_matrix.loc[pair[0], pair[1]] = 1
        adjacency_matrix.loc[pair[1], pair[0]] = 1

    n_comp, assignments = connected_components(adjacency_matrix.values)
    final_components = []
    for i in range(n_comp):
        comp = list(np.array(pair_order)[assignments == i])
        final_components.append(comp)

    final_clusters = [{seq: sequence_dict[seq]} for seq in no_match]

    for comp in final_components:
        new_seq_cluster = dict()
        for i in comp:
            new_seq_cluster[sequences[i]] = sequence_dict[sequences[i]]
        final_clusters.append(new_seq_cluster)

    return final_clusters


def merge_cluster_counts(seq_clusters, flanktrie):

    seq_clusters = [list(cluster.keys()) for cluster in seq_clusters]
    out = []
    for clust in seq_clusters:

        subtrie = flanktrie.make_subtrie(clust)
        total_words = subtrie.total_words
        out.append(total_words)

    return out


def get_cluster_consensus_seqs(seq_clusters, min_count_consensus):
    consensus_seqs = list()

    for clust in seq_clusters:
        seqs = [seq for seq in clust]
        quals = [clust[seq][0] for seq in clust]
        counts = [clust[seq][1] for seq in clust]

        mytrie = flanktrie.Trie()
        mytrie.load_words(seqs, quals, counts)
        consensus = mytrie.make_consensus_word(min_count_consensus)

        consensus_seqs.append(consensus)

    return consensus_seqs


def count_runthrough_reads(flanks, bam, min_qual, min_alignment_inner_length):
    unique_sites = flanks[['contig', 'pos']].drop_duplicates()

    runthrough_counts = []
    for index, row in unique_sites.iterrows():
        runthrough_counts.append(pysamtools.count_runthrough_reads(bam, row['contig'], row['pos'], min_qual,
                                                                   min_alignment_inner_length))

    unique_sites.loc[:,'runthrough_count'] = runthrough_counts
    outflanks = flanks.merge(unique_sites, how='left')
    return outflanks


def filter_lone_flanks(flanks, max_direct_repeat_length=30):

    column_names = ['contig']
    keep_indices = set()
    for index_5p, row in flanks.iterrows():

        if row.orient != 'R':
            continue

        contig, pos_5p = row['contig'], row['pos']

        min_pos = pos_5p - max_direct_repeat_length

        candidate_pairs = flanks.query('contig == @contig & pos > @min_pos & pos < @pos_5p & orient == "L"')

        if candidate_pairs.shape[0] > 0:
            keep_indices.add(index_5p)
            for index_3p, row in candidate_pairs.iterrows():
                keep_indices.add(index_3p)

    outpairs = flanks.ix[list(keep_indices)]
    return outpairs


def _findflanks(bamfile, min_softclip_length, min_softclip_count, min_count_consensus, min_alignment_quality,
                min_softclip_ratio, min_alignment_inner_length, output_file):
    bam = pysam.AlignmentFile(bamfile, 'rb')

    logger.info("Finding all softclipped sites...")
    softclipped_sites = get_softclipped_sites(bam, min_softclip_length, min_softclip_count, min_alignment_quality,
                                              min_alignment_inner_length, min_softclip_ratio)

    if softclipped_sites.shape[0] > 0:

        logger.info("Determining consensus flank sequences from softclipped reads...")
        flank_sequences = determine_flank_sequence(
            softclipped_sites, bam, min_softclip_length, min_softclip_count, min_count_consensus,
            min_alignment_quality, min_alignment_inner_length, min_softclip_ratio)

        # Getting the final columns
        flank_sequences = flank_sequences.loc[:,['contig', 'pos', 'orient', 'softclip_count', 'runthrough_count',
                                                'consensus_seq']]
        flank_sequences = flank_sequences.sort_values(['contig', 'pos', 'orient'])
        flank_sequences.reset_index(drop=True)
        flank_sequences.index = flank_sequences.index+1
        flank_sequences.index.name = 'flank_id'

    else:
        flank_sequences = pd.DataFrame(columns=['contig', 'pos', 'orient', 'softclip_count', 'runthrough_count', 'consensus_seq'])
        flank_sequences.index.name = 'flank_id'


    if output_file:
        logger.info("Saving results to file %s" % output_file)
        flank_sequences.to_csv(output_file, sep='\t')

    return flank_sequences


if __name__ == '__main__':
    findflanks()