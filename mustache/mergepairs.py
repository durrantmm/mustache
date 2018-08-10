import warnings
warnings.filterwarnings("ignore")
import sys
import click
import pygogo as gogo
import pandas as pd
import numpy as np
import itertools
from jellyfish import levenshtein_distance
from scipy.sparse.csgraph import connected_components

verbose=True
logger = gogo.Gogo(__name__, verbose=verbose).logger

def read_in_pairsfiles(pairsfiles):
    complete_table = None

    for f in pairsfiles:
        df = pd.read_csv(f, sep='\t')
        df.insert(loc=0, column='file', value=[f] * df.shape[0])
        if complete_table is None:
            complete_table = df
        else:
            complete_table = pd.concat([complete_table, df])

    return complete_table


def get_sequence_clusters(sequences, perc_similarity=0.9):

    if len(sequences) == 1:
        return [[0]]

    no_match = set(range(len(sequences)))
    has_match = set()

    pairwise_matches = set()
    for i, j in itertools.combinations(range(len(sequences)), r=2):

        seq1, seq2 = sequences[i], sequences [j]
        minlength = min([len(seq1), len(seq2)])
        dist = levenshtein_distance(seq1[:minlength], seq2[:minlength])
        similarity = 1 - (dist / minlength)

        if similarity >= perc_similarity:
            has_match.add(i)
            has_match.add(j)
            if i in no_match: no_match.remove(i)
            if j in no_match: no_match.remove(j)
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


    final_clusters = [[i] for i in no_match] + final_components

    return final_clusters


def _mergepairs(pairsfiles, output_file=None, similarity=0.9, comparison_sequence_length=50):

    pairs = read_in_pairsfiles(pairsfiles)
    pairs = pairs[['contig', 'pos_5p', 'pos_3p', 'seq_5p', 'seq_3p', 'blast_IS_family_5p', 'blast_IS_family_3p']].drop_duplicates()

    final_pairs = None
    sites = pairs[['contig', 'pos_5p', 'pos_3p']].drop_duplicates()
    logger.info("Beginning MERGEPAIRS on %d pairs across %d unique sites..." % (pairs.shape[0], sites.shape[0]))
    count = 0
    for index, row in sites.iterrows():
        count += 1
        if count % 100 == 0:
            logger.info("Processed %d sites, %d merged insertion sequences remain..." % (count, final_pairs.shape[0]))
        contig, pos5p, pos3p = row['contig'], row['pos_5p'], row['pos_3p']
        logger.info('Merging sequences found at site %s' % '-'.join([':'.join([contig, str(pos5p)]), str(pos3p)]))
        pairs_at_site = pairs.query("contig == @contig & pos_5p == @pos5p & pos_3p == @pos3p")

        if pairs_at_site.shape[0] > 1:

            seqs5p = list(pairs_at_site['seq_5p'])
            cropped_seqs5p = [seq[:comparison_sequence_length] if len(seq) > comparison_sequence_length else seq for seq in seqs5p]
            seq5p_clusters = get_sequence_clusters(cropped_seqs5p, similarity)
            seq5p_cluster_seqs = [list(np.array(seqs5p)[clust]) for clust in seq5p_clusters]
            seq5p_cluster_dict = {seq: max(clust, key=len) for clust in seq5p_cluster_seqs for seq in clust}
            pairs_at_site.loc[:, 'seq_5p'] = pairs_at_site.apply(lambda row, seqdict=seq5p_cluster_dict: seqdict[row['seq_5p']], axis=1)

            is_families_5p = list(pairs_at_site['blast_IS_family_5p'])
            seq5p_cluster_is_families = [list(np.array(is_families_5p)[clust]) for clust in seq5p_clusters]
            isfamily5p_cluster_dict = {family: max(set(clust), key=clust.count) for clust in seq5p_cluster_is_families for family in clust}
            pairs_at_site.loc[:, 'blast_IS_family_5p'] = pairs_at_site.apply(lambda row, seqdict=isfamily5p_cluster_dict: seqdict[row['blast_IS_family_5p']], axis=1)


            seqs3p = [seq[::-1] for seq in list(pairs_at_site['seq_3p'])]
            cropped_seqs3p = [seq[:comparison_sequence_length] if len(seq) > comparison_sequence_length else seq for seq in seqs3p]
            seq3p_clusters = get_sequence_clusters(cropped_seqs3p, similarity)
            seq3p_cluster_seqs = [list(np.array(seqs3p)[clust]) for clust in seq3p_clusters]
            seq3p_cluster_dict = {seq[::-1]: max(clust, key=len)[::-1] for clust in seq3p_cluster_seqs for seq in clust}
            pairs_at_site.loc[:, 'seq_3p'] = pairs_at_site.apply(lambda row, seqdict=seq3p_cluster_dict: seqdict[row['seq_3p']], axis=1)

            is_families_3p = list(pairs_at_site['blast_IS_family_3p'])
            seq3p_cluster_is_families = [list(np.array(is_families_3p)[clust]) for clust in seq3p_clusters]
            isfamily3p_cluster_dict = {family: max(set(clust), key=clust.count) for clust in seq3p_cluster_is_families for family in clust}
            pairs_at_site.loc[:, 'blast_IS_family_3p'] = pairs_at_site.apply(lambda row, seqdict=isfamily3p_cluster_dict: seqdict[row['blast_IS_family_3p']], axis=1)


            pairs_at_site.drop_duplicates(inplace=True)
            if pairs_at_site.shape[0] > 1:
                if len(set(list(pairs_at_site['seq_5p']))) == 1:
                    longest_seq3p = max(list(pairs_at_site['seq_3p']), key=len)
                    pairs_at_site.loc[:, 'seq_3p'] = [longest_seq3p]*pairs_at_site.shape[0]
                if len(set(list(pairs_at_site['seq_3p']))) == 1:
                    longest_seq5p = max(list(pairs_at_site['seq_5p']), key=len)
                    pairs_at_site.loc[:, 'seq_5p'] = [longest_seq5p]*pairs_at_site.shape[0]
                pairs_at_site.drop_duplicates(inplace=True)

        if final_pairs is None:
            final_pairs = pairs_at_site
        else:
            final_pairs = pd.concat([final_pairs, pairs_at_site])

    logger.info("After merging, %d unique sites remain" % final_pairs.shape[0])
    if output_file:
        logger.info("Saving results to file %s" % output_file)
        final_pairs.to_csv(output_file, sep='\t', index=False)

    return final_pairs





@click.command()
@click.argument('pairsfiles', nargs=-1, type=click.Path(exists=True))
@click.option('--output_file', '-o', default='mustache.mergepairs.tsv', help="The output file to save the results.")
@click.option('--similarity', '-s', default=0.8, help="The similarity cutoff used to merge sites.")
@click.option('--comparison_sequence_length', '-s', default=50, help="The similarity cutoff used to merge sites.")
def mergepairs(pairsfiles, output_file, similarity, comparison_sequence_length):
    if len(pairsfiles) == 0:
        logger.info("Please specify 1 or more pairs files to merge...")
        sys.exit()
    _mergepairs(pairsfiles, output_file, similarity, comparison_sequence_length)


if __name__ == '__main__':
    mergepairs()