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
verbose=True
logger = gogo.Gogo(__name__, verbose=verbose).logger


def get_softclipped_sites(bam_file, min_softclip_length=4, min_softclip_count=10, min_alignment_quality=20):

    column_names = ["contig", "pos", "orient", "count", "meets_minlength"]
    column_index_key = dict(list(zip(column_names, range(len(column_names)))))

    soft_clips = dict()

    read_count = 0
    for read in bam_file.fetch():
        read_count += 1
        if read_count % 100000 == 0:
            logger.info("\tAfter checking %d reads, %d softclipped sites found..." % (read_count, len(soft_clips)))

        if read.mapping_quality < min_alignment_quality:
            continue

        if sctools.is_right_softclipped_strict(read):

            meets_minlength = False
            if sctools.right_softclip_length(read) >= min_softclip_length:
                meets_minlength = True

            contig, pos = sctools.right_softclipped_site_lenient(read)
            orient = 'R'
            index = (contig, pos, orient)

            if index in soft_clips:
                soft_clips[index][column_index_key['count']] += 1
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
                soft_clips[index][column_index_key['count']] += 1
                if soft_clips[index][column_index_key['meets_minlength']] == False and meets_minlength:
                    soft_clips[index][column_index_key['meets_minlength']] = meets_minlength
            else:
                new_site = [contig, pos, orient, 1, meets_minlength]
                soft_clips[index] = new_site

    soft_clips = pd.DataFrame.from_dict(soft_clips, orient='index', columns=column_names).sort_index()
    soft_clips.query('meets_minlength == True & count >= @min_softclip_count', inplace=True)

    logger.info("After filtering by minimum softclip length of %d, "
                "and a minimum softclip count of %d, %d sites remaining..." % (min_softclip_length, min_softclip_count, soft_clips.shape[0]))

    return soft_clips


def determine_flank_sequence(softclipped_sites, bam_file, min_softclip_length=4, min_softclip_count=10, min_alignment_quality=20):

    column_names = ["contig", "pos", "orient", "count", "consensus_seq_length", "consensus_seq"]

    final_flank_seqs = dict()
    site_count = 0
    for index, row in softclipped_sites.iterrows():
        site_count += 1
        if site_count % 1000 == 0:
            logger.info("\tAfter processing %d sites, %d unique softclipped flanks found..." % (site_count, len(final_flank_seqs)))

        contig, pos, orient = index


        softclipped_seqs, softclipped_qualities = pysamtools.get_softclipped_seqs_qualities(bam_file, contig, pos,
                                                                                            orient, min_alignment_quality)

        flanktrie = get_flank_sequence_trie(softclipped_seqs, softclipped_qualities)
        seq_dict = dict(flanktrie.traverse_both())

        seq_clusters = get_sequence_clusters(seq_dict)

        consensus_seqs = get_cluster_consensus_seqs(seq_clusters)

        if orient == 'L':
            consensus_seqs = [seq[::-1] for seq in consensus_seqs]

        merged_counts = merge_cluster_counts(seq_clusters, flanktrie)


        for seq, count in zip(consensus_seqs, merged_counts):
            final_flank_seqs[len(final_flank_seqs)+1] = [contig, pos, orient, count, len(seq), seq]

    final_flank_seqs = pd.DataFrame.from_dict(final_flank_seqs, orient='index', columns=column_names).sort_index()
    final_flank_seqs.query('consensus_seq_length >= @min_softclip_length & count >= @min_softclip_count', inplace=True)

    logger.info("After filtering by minimum softclip length of %d, "
                "and a minimum softclip count of %d, %d sites remaining..." % (
                min_softclip_length, min_softclip_count, final_flank_seqs.shape[0]))
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


def get_cluster_consensus_seqs(seq_clusters):
    consensus_seqs = list()

    for clust in seq_clusters:
        seqs = [seq for seq in clust]
        quals = [clust[seq] for seq in clust]

        max_len = max(map(len, seqs))

        consensus = ''
        for i in range(max_len):
            bestbase = None
            bestqual = -1

            for j in range(len(seqs)):
                word = seqs[j]
                qual = quals[j]

                if i >= len(word):
                    continue

                base = word[i]
                q = qual[i]

                if q > bestqual:
                    bestqual = q
                    bestbase = base

            consensus += bestbase


        consensus_seqs.append(consensus)

    return consensus_seqs

def _findflanks(bamfile, min_softclip_length, min_softclip_count, min_alignment_quality, output_file):
    bam = pysam.AlignmentFile(bamfile, 'rb')

    logger.info("Finding all softclipped sites...")
    softclipped_sites = get_softclipped_sites(bam, min_softclip_length, min_softclip_count, min_alignment_quality)
    logger.info("Determining consensus flank sequences from softclipped reads...")
    flank_sequences = determine_flank_sequence(softclipped_sites, bam, min_softclip_length, min_softclip_count,
                                               min_alignment_quality)

    if output_file:
        logger.info("Saving results to file %s" % output_file)
        flank_sequences.to_csv(output_file, sep='\t', index=False)

    return flank_sequences

@click.command()
@click.argument('bamfile', type=click.Path(exists=True))
@click.option('--output_file', '-o', default='mustache.findflanks.tsv', help="The output file to save the results.")
@click.option('--min_softclip_length', '-minlen', default=4, help="For a softclipped site to be considered, there must be at least one softclipped read of this length.")
@click.option('--min_softclip_count', '-mincount', default=10, help="For a softclipped site to be considered, there must be at least this many softclipped reads at the site.")
@click.option('--min_alignment_quality', '-minq', default=20, help="For a read to be considered, it must meet this alignment quality cutoff.")
def findflanks(bamfile, min_softclip_length, min_softclip_count, min_alignment_quality, output_file=None):
    _findflanks(bamfile, min_softclip_length, min_softclip_count, min_alignment_quality, output_file)



if __name__ == '__main__':
    findflanks()
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