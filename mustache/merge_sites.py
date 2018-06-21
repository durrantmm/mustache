import sys
from collections import defaultdict, OrderedDict
import pandas as pd
from mustache import output, misc
import numpy as np
import click
from mustache.assemble_tools import MinimusAssembler, Minimus2Merger
from mustache.alignment_tools import merge_flank_assemblies
from os.path import join
from tqdm import tqdm
tqdm.monitor_interval = 0
import itertools
from scipy.sparse.csgraph import connected_components
from Bio.pairwise2 import format_alignment
from Bio import pairwise2

def merge(insertseq_files, filt_min_merged_length, filt_min_flank_length, pairwise_ident_cutoff,
          min_softclip_pair_distance, max_softclip_pair_distance, output_prefix, lenient, outdir):

    complete_table = None

    for f in insertseq_files:
        df = pd.read_csv(f, sep='\t')
        df.insert(loc=0, column='file', value=[f]*df.shape[0])
        if complete_table is None:
            complete_table = df
        else:
            complete_table = pd.concat([complete_table, df])

    click.echo("%d insertion sequences to process...\n" % complete_table.shape[0])

    merged_right_sites = merge_right_sites(complete_table, filt_min_flank_length, pairwise_ident_cutoff)
    click.echo("After processing right sites, %d merged sequences remain...\n" % merged_right_sites.shape[0])

    merged_left_sites = merge_left_sites(complete_table, filt_min_flank_length, pairwise_ident_cutoff)
    click.echo("After processing left sites, %d merged sequences remain...\n" % merged_left_sites.shape[0])



    left_right_sites, merged_left_right_sites = merge_left_right_sites(complete_table,
                                                                       merged_right_sites, merged_left_sites,
                                                                       min_softclip_pair_distance,
                                                                       max_softclip_pair_distance, outdir)

    all_merged_sorted = complete_table.query('orientation == "M"').sort_values(['contig', 'left_site', 'right_site'])
    all_merged_sorted = all_merged_sorted[['contig', 'left_site', 'right_site', 'orientation', 'partner_site', 'assembly_length', 'assembly']]
    all_merged_sorted = pd.concat([all_merged_sorted, merged_left_right_sites]).sort_values(['contig', 'site']).reset_index(drop=True)

    final_merged_full_sites = merge_full_sites_strict(all_merged_sorted, filt_min_merged_length, pairwise_ident_cutoff)
    if lenient:
        final_merged_full_sites = merge_full_sites_lenient(final_merged_full_sites, filt_min_merged_length)

    click.echo("After processing full sites, %d merged sequences remain...\n" % final_merged_full_sites.shape[0])

    final_left_right_sites = compare_flanks_to_merged_sites(left_right_sites, final_merged_full_sites, pairwise_ident_cutoff)

    final_output_table = pd.concat([final_left_right_sites, final_merged_full_sites]).sort_values(['contig', 'site']).reset_index(drop=True)[
        ['contig', 'left_site', 'right_site', 'orientation', 'partner_site', 'assembly_length', 'assembly']].drop_duplicates()

    click.echo("After processing %d initial input sites, %d final sites remain...\n" %
               (complete_table.shape[0], final_output_table.shape[0]))

    click.echo("Writing final results to output table and FASTA file...\n")
    output.write_final_dataframe(final_output_table, outdir, output_prefix, '.merged_insertion_seqs.tsv')
    output.write_final_merged_dataframe_to_fasta(final_output_table, outdir, output_prefix)
    click.echo("Merge complete...")

    sys.exit()


def merge_full_sites_strict(full_sorted, filt_min_merged_length, pairwise_ident_cutoff):

    click.echo("Filtering out all merged sequences shorter than %d base pairs..." % filt_min_merged_length)
    full_sorted.query("assembly_length >= @filt_min_merged_length", inplace=True)

    unique_sites = full_sorted[['contig', 'left_site', 'right_site']].drop_duplicates()
    click.echo('Processing %d complete insertion sites...' % unique_sites.shape[0])

    with tqdm(total=unique_sites.shape[0], desc='MERGE FULL SITES STRICT') as bar:

        out_dict = OrderedDict([
            ('contig', []),
            ('site', []),
            ('left_site', []),
            ('right_site', []),
            ('orientation', []),
            ('partner_site', []),
            ('assembly_length', []),
            ('assembly', [])
        ])

        for index, row in unique_sites.iterrows():
            bar.update(1)
            contig, left_site, right_site = row['contig'], row['left_site'], row['right_site']
            matched_sites = full_sorted.query('contig == @contig & left_site == @left_site & right_site == @right_site')

            unique_assemblies = np.array(list(set(list(matched_sites['assembly']))))

            if len(unique_assemblies) == 1:
                out_seq = list(unique_assemblies)[0]

                out_dict['contig'].append(contig)
                out_dict['site'].append(left_site)
                out_dict['left_site'].append(left_site)
                out_dict['right_site'].append(right_site)
                out_dict['orientation'].append('M')
                out_dict['partner_site'].append(-1)
                out_dict['assembly_length'].append(len(out_seq))
                out_dict['assembly'].append(out_seq)

            else:
                same_seq_length_sets = get_same_sequence_length_sets(unique_assemblies)

                for unique_assemblies in same_seq_length_sets:

                    out_sequences = np.array(list(unique_assemblies))
                    out_sequences_prev_length = -1

                    while len(out_sequences) != out_sequences_prev_length:

                        has_match = set()
                        pairwise_matches = set()
                        for i, j in itertools.combinations(range(len(out_sequences)), r=2):
                            pairwise_ident = get_perc_identity(out_sequences[i], out_sequences[j])
                            if pairwise_ident > pairwise_ident_cutoff:
                                has_match.add(i)
                                has_match.add(j)
                                pairwise_matches.add((i, j))

                        new_out_sequences = [out_sequences[i] for i in range(len(out_sequences)) if i not in has_match]

                        pairwise_clusters = get_pairwise_clusters(pairwise_matches, has_match)

                        for cluster in pairwise_clusters:
                            cluster_seqs = out_sequences[cluster]
                            longest_seq = max(cluster_seqs, key=len)
                            new_out_sequences.append(longest_seq)

                        out_sequences_prev_length = len(out_sequences)
                        out_sequences = np.array(list(new_out_sequences))

                    for seq in out_sequences:
                        out_dict['contig'].append(contig)
                        out_dict['site'].append(left_site)
                        out_dict['left_site'].append(left_site)
                        out_dict['right_site'].append(right_site)
                        out_dict['orientation'].append('M')
                        out_dict['partner_site'].append(-1)
                        out_dict['assembly_length'].append(len(seq))
                        out_dict['assembly'].append(seq)

    merged_full_sites = pd.DataFrame(out_dict).drop_duplicates()
    merged_full_sites.reset_index(drop=True, inplace=True)
    return merged_full_sites


def merge_full_sites_lenient(full_sorted, filt_min_merged_length, pairwise_indel_cutoff=0.15):

    click.echo("Filtering out all merged sequences shorter than %d base pairs..." % filt_min_merged_length)
    full_sorted.query("assembly_length >= @filt_min_merged_length", inplace=True)

    unique_sites = full_sorted[['contig', 'left_site', 'right_site']].drop_duplicates()
    click.echo('Processing %d complete insertion sites...' % unique_sites.shape[0])

    with tqdm(total=unique_sites.shape[0], desc='MERGE FULL SITES LENIENT') as bar:

        out_dict = OrderedDict([
            ('contig', []),
            ('site', []),
            ('left_site', []),
            ('right_site', []),
            ('orientation', []),
            ('partner_site', []),
            ('assembly_length', []),
            ('assembly', [])
        ])

        for index, row in unique_sites.iterrows():
            bar.update(1)
            contig, left_site, right_site = row['contig'], row['left_site'], row['right_site']
            matched_sites = full_sorted.query('contig == @contig & left_site == @left_site & right_site == @right_site')

            unique_assemblies = np.array(list(set(list(matched_sites['assembly']))))

            if len(unique_assemblies) == 1:
                out_seq = list(unique_assemblies)[0]

                out_dict['contig'].append(contig)
                out_dict['site'].append(left_site)
                out_dict['left_site'].append(left_site)
                out_dict['right_site'].append(right_site)
                out_dict['orientation'].append('M')
                out_dict['partner_site'].append(-1)
                out_dict['assembly_length'].append(len(out_seq))
                out_dict['assembly'].append(out_seq)

            else:

                out_sequences = np.array(list(unique_assemblies))
                out_sequences_prev_length = -1

                while len(out_sequences) != out_sequences_prev_length:

                    has_match = set()
                    pairwise_matches = set()
                    for i, j in itertools.combinations(range(len(out_sequences)), r=2):
                        pairwise_indel = get_perc_indels_global(out_sequences[i], out_sequences[j])
                        if pairwise_indel < pairwise_indel_cutoff:
                            has_match.add(i)
                            has_match.add(j)
                            pairwise_matches.add((i, j))

                    new_out_sequences = [out_sequences[i] for i in range(len(out_sequences)) if i not in has_match]

                    pairwise_clusters = get_pairwise_clusters(pairwise_matches, has_match)

                    for cluster in pairwise_clusters:

                        cluster_seqs = out_sequences[cluster]
                        consensus_sequence = get_consensus_sequence(cluster_seqs)
                        new_out_sequences.append(consensus_sequence)

                    out_sequences_prev_length = len(out_sequences)
                    out_sequences = np.array(list(new_out_sequences))

                for seq in out_sequences:
                    out_dict['contig'].append(contig)
                    out_dict['site'].append(left_site)
                    out_dict['left_site'].append(left_site)
                    out_dict['right_site'].append(right_site)
                    out_dict['orientation'].append('M')
                    out_dict['partner_site'].append(-1)
                    out_dict['assembly_length'].append(len(seq))
                    out_dict['assembly'].append(seq)

    merged_full_sites = pd.DataFrame(out_dict).drop_duplicates()
    merged_full_sites.reset_index(drop=True, inplace=True)
    return merged_full_sites

def merge_right_sites(complete_table, filt_min_flank_length, pairwise_ident_cutoff):
    right_sorted = complete_table.query('orientation == "R"').sort_values(['contig', 'right_site'], ascending=True)
    right_sorted.reset_index(drop=True, inplace=True)

    click.echo("Filtering out all right flank sequences shorter than %d base pairs..." % filt_min_flank_length)
    right_sorted.query("assembly_length >= @filt_min_flank_length", inplace=True)

    unique_sites = pd.DataFrame(right_sorted[['contig', 'right_site']].drop_duplicates())
    click.echo('Processing %d right flank insertion sites...' % unique_sites.shape[0])


    with tqdm(total=unique_sites.shape[0], desc='MERGE RIGHT SITES') as bar:

        out_dict = OrderedDict([
            ('contig', []),
            ('site', []),
            ('left_site', []),
            ('right_site', []),
            ('orientation', []),
            ('assembly_length', []),
            ('assembly', [])
        ])

        for index, row in unique_sites.iterrows():
            bar.update(1)
            contig, right_site = row['contig'], row['right_site']
            matched_sites = right_sorted.query('contig == @contig & right_site == @right_site')

            unique_assemblies = np.array(list(set(list(matched_sites['assembly']))))

            if len(unique_assemblies) == 1:
                out_seq = list(unique_assemblies)[0]

                out_dict['contig'].append(contig)
                out_dict['site'].append(right_site)
                out_dict['left_site'].append(-1)
                out_dict['right_site'].append(right_site)
                out_dict['orientation'].append('R')
                out_dict['assembly_length'].append(len(out_seq))
                out_dict['assembly'].append(out_seq)

            else:

                out_sequences = np.array(list(unique_assemblies))
                out_sequences_prev_length = -1

                while len(out_sequences) != out_sequences_prev_length:

                    has_match = set()
                    pairwise_matches = set()
                    for i, j in itertools.combinations(range(len(out_sequences)), r=2):
                        pairwise_ident = get_perc_identity(out_sequences[i], out_sequences[j])
                        if pairwise_ident > pairwise_ident_cutoff:
                            has_match.add(i)
                            has_match.add(j)
                            pairwise_matches.add((i, j))

                    new_out_sequences = [out_sequences[i] for i in range(len(out_sequences)) if i not in has_match]

                    pairwise_clusters = get_pairwise_clusters(pairwise_matches, has_match)

                    for cluster in pairwise_clusters:
                        cluster_seqs = out_sequences[cluster]
                        consensus_sequence = get_consensus_sequence(cluster_seqs)
                        new_out_sequences.append(consensus_sequence)

                    out_sequences_prev_length = len(out_sequences)
                    out_sequences = np.array(list(new_out_sequences))

                for seq in out_sequences:
                    out_dict['contig'].append(contig)
                    out_dict['site'].append(right_site)
                    out_dict['left_site'].append(-1)
                    out_dict['right_site'].append(right_site)
                    out_dict['orientation'].append('R')
                    out_dict['assembly_length'].append(len(seq))
                    out_dict['assembly'].append(seq)

    merged_right_sites = pd.DataFrame(out_dict)
    merged_right_sites.reset_index(drop=True, inplace=True)
    return merged_right_sites


def merge_left_sites(complete_table, filt_min_flank_length, pairwise_ident_cutoff):
    left_sorted = complete_table.query('orientation == "L"').sort_values(['contig', 'left_site'], ascending=True)
    left_sorted.reset_index(drop=True, inplace=True)

    click.echo("Filtering out all right flank sequences shorter than %d base pairs..." % filt_min_flank_length)
    left_sorted.query("assembly_length >= @filt_min_flank_length", inplace=True)

    unique_sites = pd.DataFrame(left_sorted[['contig', 'left_site']].drop_duplicates())
    click.echo('Processing %d left flank insertion sites...' % unique_sites.shape[0])

    with tqdm(total=unique_sites.shape[0], desc='MERGE LEFT SITES') as bar:
        out_dict = OrderedDict([
            ('contig', []),
            ('site', []),
            ('left_site', []),
            ('right_site', []),
            ('orientation', []),
            ('assembly_length', []),
            ('assembly', [])
        ])


        for index, row in unique_sites.iterrows():
            bar.update(1)
            contig, left_site = row['contig'], row['left_site']
            matched_sites = left_sorted.query('contig == @contig & left_site == @left_site')

            unique_assemblies = np.array(list(set(list(matched_sites['assembly']))))

            if len(unique_assemblies) == 1:
                out_seq = list(unique_assemblies)[0]

                out_dict['contig'].append(contig)
                out_dict['site'].append(left_site)
                out_dict['left_site'].append(left_site)
                out_dict['right_site'].append(-1)
                out_dict['orientation'].append('L')
                out_dict['assembly_length'].append(len(out_seq))
                out_dict['assembly'].append(out_seq)

            else:

                out_sequences = np.array(list(unique_assemblies))
                out_sequences_prev_length = -1

                while len(out_sequences) != out_sequences_prev_length:

                    has_match = set()
                    pairwise_matches = set()
                    for i, j in itertools.combinations(range(len(out_sequences)), r=2):
                        pairwise_ident = get_perc_identity(out_sequences[i], out_sequences[j], is_reversed=True)
                        if pairwise_ident > pairwise_ident_cutoff:
                            has_match.add(i)
                            has_match.add(j)
                            pairwise_matches.add((i, j))

                    new_out_sequences = [out_sequences[i] for i in range(len(out_sequences)) if i not in has_match]

                    pairwise_clusters = get_pairwise_clusters(pairwise_matches, has_match)

                    for cluster in pairwise_clusters:
                        cluster_seqs = out_sequences[cluster]
                        consensus_sequence = get_consensus_sequence(cluster_seqs, is_reversed=True)
                        new_out_sequences.append(consensus_sequence)

                    out_sequences_prev_length = len(out_sequences)
                    out_sequences = np.array(list(new_out_sequences))

                for seq in out_sequences:
                    out_dict['contig'].append(contig)
                    out_dict['site'].append(left_site)
                    out_dict['left_site'].append(left_site)
                    out_dict['right_site'].append(-1)
                    out_dict['orientation'].append('L')
                    out_dict['assembly_length'].append(len(seq))
                    out_dict['assembly'].append(seq)

    merged_left_sites = pd.DataFrame(out_dict)
    merged_left_sites.reset_index(drop=True, inplace=True)

    return merged_left_sites


def merge_left_right_sites(complete_table, merged_right_sites, merged_left_sites, min_softclip_pair_distance, max_softclip_pair_distance, outdir):

    left_site_pairs = complete_table.query('left_site != -1 & orientation != "M"')[['contig', 'left_site', 'partner_site']].drop_duplicates()
    right_site_pairs = complete_table.query('right_site != -1 & orientation != "M"')[['contig', 'right_site', 'partner_site']].drop_duplicates()

    all_flank_sites = pd.concat([merged_right_sites, merged_left_sites]).sort_values(['contig', 'site']).reset_index(drop=True)

    left_sites = all_flank_sites.query('left_site > -1')
    right_sites = all_flank_sites.query('right_site > -1')

    with tqdm(total=2, desc='PAIRING FLANKS') as bar:
        left_sites = pd.merge(left_sites, left_site_pairs, how='inner')
        bar.update(1)
        right_sites = pd.merge(right_sites, right_site_pairs, how='inner')
        bar.update(1)


    left_right_sites = pd.concat([left_sites, right_sites]).sort_values(['contig', 'site']).reset_index(drop=True).drop_duplicates()
    left_sites = left_right_sites.query('left_site > -1')

    with tqdm(total=left_sites.shape[0], desc='ATTEMPTING FLANK MERGE') as bar:

        out_dict = OrderedDict([ ('contig', []), ('site', []), ('left_site', []), ('right_site', []),
            ('orientation', []), ('partner_site', []), ('assembly_length', []), ('assembly', [])])

        out_dict_merged = OrderedDict([('contig', []), ('site', []), ('left_site', []), ('right_site', []),
                                ('orientation', []), ('partner_site', []), ('assembly_length', []), ('assembly', [])])

        for index1, row1 in left_sites.iterrows():
            bar.update(1)
            left_partner_site = row1['partner_site']
            left_assembly = row1['assembly']

            right_partners = left_right_sites.query('contig == "{contig}" & '
                                                    'right_site == @left_partner_site'.format(contig=row1['contig']))

            for index2, row2 in right_partners.iterrows():
                right_assembly = row2['assembly']
                merged_assembly = merge_flank_assemblies(right_assembly, left_assembly, verbose=False)

                if merged_assembly is None:

                    out_dict['contig'].append(row1['contig'])
                    out_dict['contig'].append(row2['contig'])
                    out_dict['site'].append(row1['site'])
                    out_dict['site'].append(row2['site'])
                    out_dict['left_site'].append(row1['left_site'])
                    out_dict['left_site'].append(row2['left_site'])
                    out_dict['right_site'].append(row1['right_site'])
                    out_dict['right_site'].append(row2['right_site'])
                    out_dict['orientation'].append(row1['orientation'])
                    out_dict['orientation'].append(row2['orientation'])
                    out_dict['partner_site'].append(row1['partner_site'])
                    out_dict['partner_site'].append(row2['partner_site'])
                    out_dict['assembly_length'].append(row1['assembly_length'])
                    out_dict['assembly_length'].append(row2['assembly_length'])
                    out_dict['assembly'].append(row1['assembly'])
                    out_dict['assembly'].append(row2['assembly'])

                else:
                    out_dict_merged['contig'].append(row1['contig'])
                    out_dict_merged['site'].append(row1['site'])
                    out_dict_merged['left_site'].append(row1['left_site'])
                    out_dict_merged['right_site'].append(row2['right_site'])
                    out_dict_merged['orientation'].append("M")
                    out_dict_merged['partner_site'].append(-1)
                    out_dict_merged['assembly_length'].append(len(merged_assembly))
                    out_dict_merged['assembly'].append(merged_assembly)

        final_left_right_sites = pd.DataFrame(out_dict).reset_index(drop=True).sort_values(['contig', 'site']).drop_duplicates()
        final_merged_sites = pd.DataFrame(out_dict_merged).reset_index(drop=True).sort_values(['contig', 'site']).drop_duplicates()

        return final_left_right_sites.drop_duplicates(), final_merged_sites.drop_duplicates()


def compare_flanks_to_merged_sites(left_right_sites, all_merged_sorted, pairwise_ident_cutoff):

    left_sites = left_right_sites.query('left_site > -1')

    final_left_to_right_sites = None
    with tqdm(total=left_sites.shape[0], desc='COMPARING FLANKS TO MERGED SITES') as bar:
        for index1, row1 in left_sites.iterrows():

            bar.update(1)

            left_site = row1['left_site']
            left_partner_site = row1['partner_site']
            right_partners = left_right_sites.query('contig == "{contig}" & '
                                                    'right_site == @left_partner_site'.format(contig=row1['contig']))

            for index2, row2 in right_partners.iterrows():
                right_site = row2['right_site']
                merged_matches = all_merged_sorted.query('contig == "{contig}" & '
                                                         'left_site == @left_site &'
                                                         'right_site == @right_site'.format(contig=row1['contig']))


                if merged_matches.shape[0] == 0:
                    out_pair = pd.DataFrame.transpose(pd.concat([pd.DataFrame(row1), pd.DataFrame(row2)], axis=1))
                    if final_left_to_right_sites is None:
                        final_left_to_right_sites = out_pair
                    else:
                        final_left_to_right_sites = pd.concat([final_left_to_right_sites, out_pair])
                else:

                    left_assembly = row1['assembly']
                    right_assembly = row2['assembly']

                    for index3, row3 in merged_matches.iterrows():
                        merged_assembly = row3['assembly']
                        right_perc_identity = get_perc_identity(right_assembly, merged_assembly)
                        left_perc_identity = get_perc_identity(left_assembly, merged_assembly, is_reversed=True)

                        if right_perc_identity > pairwise_ident_cutoff and left_perc_identity > pairwise_ident_cutoff:
                            break
                    else:

                        out_pair = pd.DataFrame.transpose(pd.concat([pd.DataFrame(row1), pd.DataFrame(row2)], axis=1))

                        if final_left_to_right_sites is None:
                            final_left_to_right_sites = out_pair
                        else:
                            final_left_to_right_sites = pd.concat([final_left_to_right_sites, out_pair])

    return final_left_to_right_sites


def get_perc_identity(seq1, seq2, is_reversed=False):
    total = min(map(len, [seq1, seq2]))
    match = 0

    if not is_reversed:
        for i in range(total):

            if seq1[i] == seq2[i]:
                match += 1
            elif seq1[i].upper() == 'N' or seq2[i].upper() == 'N':
                match += 1
    else:
        revseq1, revseq2 = seq1[::-1], seq2[::-1]
        for i in range(total):
            if revseq1[i] == revseq2[i]:
                match += 1
            elif revseq1[i].upper() == 'N' or revseq2[i].upper() == 'N':
                match += 1

    return match/total


def get_perc_indels_global(seq1, seq2, is_reversed=False):
    alignments = pairwise2.align.globalxx(seq1, seq2)
    aln1, aln2 = alignments[0][0], alignments[0][1]


    inside_indel = False
    total_chars = 0
    total_indels = 0

    for i in range(len(aln1)):
        if aln1[i] == '-' or aln2[i] == '-':

            if not inside_indel:
                inside_indel = True
                total_chars += 1
                total_indels += 1

        elif aln1[i] != '-' and aln2[i] != '-':
            total_chars += 1
            if inside_indel:
                inside_indel = False

    perc_indels = total_indels / total_chars

    #if perc_indels > 0:
    #    print(total_indels)
    #    print(total_chars)
    #    print(perc_indels)
    #    print(format_alignment(*alignments[0]))
    return perc_indels


def get_pairwise_clusters(pairwise_matches, has_match):
    pair_order = sorted(list(has_match))
    adjacency_matrix = pd.DataFrame(np.zeros((len(pair_order), len(pair_order))), index=pair_order, columns=pair_order)

    for pair in pairwise_matches:
        adjacency_matrix.loc[pair[0], pair[1]] = 1
        adjacency_matrix.loc[pair[1], pair[0]] = 1

    n_comp, assignments = connected_components(adjacency_matrix)

    final_components = []
    for i in range(n_comp):
        comp = list(np.array(pair_order)[assignments == i])
        final_components.append(comp)

    return final_components

def get_consensus_sequence(sequences, is_reversed=False):

    max_seq_len = max(map(len, sequences))
    count_matrix = np.zeros((4, max_seq_len))
    for seq in sequences:
        if is_reversed:
            seq = seq[::-1]
        for i in range(len(seq)):
            try:
                nuci = nuc_to_num(seq[i])
                if nuci != -1:
                    count_matrix[nuci, i] += 1
            except IndexError:
                continue
    maxima = np.apply_along_axis(lambda x: num_to_nuc(np.where(x == max(x))[0]) if sum(x == max(x)) == 1 else 'N',
                                 axis=0, arr=count_matrix)

    consensus = ''.join(list(maxima))

    if is_reversed:
        consensus = consensus[::-1]
    return consensus

def nuc_to_num(nuc):
    if nuc.upper() == 'A':
        return 0
    elif nuc.upper() == 'C':
        return 1
    elif nuc.upper() == 'G':
        return 2
    elif nuc.upper() == 'T':
        return 3
    else:
        return -1

def num_to_nuc(nuc):
    if nuc == 0:
        return 'A'
    elif nuc == 1:
        return 'C'
    elif nuc == 2:
        return 'G'
    elif nuc == 3:
        return 'T'
    else:
        print(nuc)
        sys.exit()

def get_same_sequence_length_sets(sequences):
    same_length = defaultdict(list)
    for seq in sequences:
        same_length[len(seq)].append(seq)
    sets = []
    for length in same_length:
        sets.append(same_length[length])
    return sets
