from collections import defaultdict, OrderedDict
from scipy.stats import poisson
import pandas as pd
import sys
import time
from statsmodels.sandbox.stats.multicomp import fdrcorrection0 as fdr_adjust

def run(bam_file, contig, start, stop, min_softclip_length, min_softclip_count,
        min_softclip_pair_distance, max_softclip_pair_distance, max_softclip_count_ratio_deviation):
    softclipped_sites = get_softclipped_sites(bam_file, contig, start, stop, min_softclip_length, min_softclip_count)
    softclipped_pairs = get_softclipped_pairs(softclipped_sites, min_softclip_pair_distance, max_softclip_pair_distance,
                                              max_softclip_count_ratio_deviation)
    return softclipped_pairs

def get_softclipped_sites(bam_file, contig, start, stop, min_softclip_length, min_softclip_count):
    soft_clips = defaultdict(lambda: defaultdict(set))

    for read in bam_file.fetch(contig, start, stop):

        if is_right_softclipped(read) and right_softclip_length(read) >= min_softclip_length:

            soft_clips[right_softclip_site(read)]['R'].add(read.query_name)

        if is_left_softclipped(read) and left_softclip_length(read) >= min_softclip_length:

            soft_clips[left_softclip_site(read)]['L'].add(read.query_name)

    soft_clip_counts = get_softclip_counts(soft_clips)
    soft_clip_counts_filtered = filter_softclip_counts(soft_clip_counts, min_softclip_count)

    return soft_clip_counts_filtered

def get_softclipped_pairs(softclipped_sites, min_softclip_pair_distance, max_softclip_pair_distance,
                          max_softclip_count_ratio_deviation):
    df = softclipped_sites
    final_pairs = OrderedDict([('contig',[]), ('left_site',[]), ('right_site',[]), ('distance', []), ('left_count',[]), ('right_count',[])])

    for contig in df.contig.unique():
        df_contig = df.query("contig == '%s'" % contig).reset_index(drop=True)

        i = 0
        while i < (len(df_contig)-1):

            L_count_p1 = df_contig.loc[i,'L_count']
            R_count_p2 = df_contig.loc[i+1, 'R_count']

            base_p1 = df_contig.loc[i, 'base']
            base_p2 = df_contig.loc[i + 1, 'base']

            distance = base_p2-base_p1
            if L_count_p1 > 0 and R_count_p2 > 0 and min_softclip_pair_distance <= distance <= max_softclip_pair_distance:
                count_ratio = float(L_count_p1) / (L_count_p1 + R_count_p2)
                if 0.5-max_softclip_count_ratio_deviation < count_ratio < 0.5+max_softclip_count_ratio_deviation:
                    final_pairs['contig'].append(contig)
                    final_pairs['left_site'].append(base_p1)
                    final_pairs['right_site'].append(base_p2)
                    final_pairs['distance'].append(base_p2-base_p1)
                    final_pairs['left_count'].append(L_count_p1)
                    final_pairs['right_count'].append(R_count_p2)
                    i += 2
                    continue

            i += 1

    out_df = pd.DataFrame(final_pairs)

    return out_df


def get_softclip_counts(soft_clips):

    soft_clip_counts = OrderedDict([('contig', []), ('base', []), ('L_count', []), ('R_count', [])])
    for site in soft_clips:
        left_softclip_count = len(soft_clips[site]['L'])
        right_softclip_count = len(soft_clips[site]['R'])
        soft_clip_counts['contig'].append(site[0])
        soft_clip_counts['base'].append(site[1])
        soft_clip_counts['L_count'].append(left_softclip_count)
        soft_clip_counts['R_count'].append(right_softclip_count)

    soft_clip_counts = pd.DataFrame(soft_clip_counts)

    return soft_clip_counts


def filter_softclip_counts(soft_clip_counts, min_softclip_count):
    df = soft_clip_counts
    df['L_count'] = [c if c >= min_softclip_count else 0 for c in df['L_count']]
    df['R_count'] = [c if c >= min_softclip_count else 0 for c in df['R_count']]

    df = df.query('L_count > 0 or R_count > 0')
    df = df.sort_values(['contig', 'base'])
    df = df.reset_index(drop=True)

    return df


def is_left_softclipped(read):
    return read.cigartuples[0][0] == 4

def is_right_softclipped(read):
    return read.cigartuples[-1][0] == 4

def left_softclip_length(read):
    return read.cigartuples[0][1]

def right_softclip_length(read):
    return read.cigartuples[-1][1]

def left_softclip_site(read):
    return (read.reference_name, read.get_reference_positions()[0] - 1)

def right_softclip_site(read):
    return (read.reference_name, read.get_reference_positions()[-1] + 1)