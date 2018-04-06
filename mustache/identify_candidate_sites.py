from collections import defaultdict, OrderedDict
import pandas as pd
from mustache import misc
from scipy.stats import poisson
import numpy as np
import click

def get_softclipped_sites(bam_file, contig, start, stop, min_softclip_length, min_softclip_count, min_softclip_pair_distance, max_softclip_pair_distance):

    soft_clips = defaultdict(lambda: defaultdict(set))

    contigs = misc.get_bam_contig_dict(bam_file)

    if contig:
        contigs= {contig: contigs[contig]}

    all_softclipped_sites = None
    for contig in contigs:
        click.echo("\tProcessing softclipped sites on contig %s..." % contig)
        for read in bam_file.fetch(contig, start, stop):

            if is_right_softclipped(read) and right_softclip_length(read) >= min_softclip_length:

                soft_clips[right_softclip_site(read)]['R'].add(read.query_name)

            if is_left_softclipped(read) and left_softclip_length(read) >= min_softclip_length:

                soft_clips[left_softclip_site(read)]['L'].add(read.query_name)

        click.echo("\tGetting softclip counts...")
        soft_clip_counts = get_softclip_counts(soft_clips)
        click.echo("\tFiltering softclipped sites...")
        soft_clip_counts_filtered = filter_softclip_counts(soft_clip_counts, min_softclip_count, contigs[contig], min_softclip_pair_distance, max_softclip_pair_distance)

        if all_softclipped_sites is None:
            all_softclipped_sites = soft_clip_counts_filtered
        else:
            all_softclipped_sites = pd.concat([all_softclipped_sites, soft_clip_counts_filtered])

    return all_softclipped_sites


def get_softclip_counts(soft_clips):

    soft_clip_counts = OrderedDict([('contig', []), ('site', []), ('L_count', []), ('R_count', [])])
    for site in soft_clips:
        left_softclip_count = len(soft_clips[site]['L'])
        right_softclip_count = len(soft_clips[site]['R'])
        soft_clip_counts['contig'].append(site[0])
        soft_clip_counts['site'].append(site[1])
        soft_clip_counts['L_count'].append(left_softclip_count)
        soft_clip_counts['R_count'].append(right_softclip_count)

    soft_clip_counts = pd.DataFrame(soft_clip_counts)
    soft_clip_counts = soft_clip_counts.sort_values(['contig', 'site']).reset_index(drop=True)

    return soft_clip_counts


def filter_softclip_counts(soft_clip_counts, min_softclip_count, contig_length, min_softclip_pair_distance, max_softclip_pair_distance,
                           flank_small_len=1000, flank_int_len=5000, flank_large_len=10000, pval_cutoff=1e-5):
    df = soft_clip_counts

    L_lambda_bg = df.L_count.sum() / contig_length
    R_lambda_bg = df.R_count.sum() / contig_length

    df.loc[:,'L_lambda_bg'] = L_lambda_bg
    df.loc[:,'R_lambda_bg'] = R_lambda_bg

    L_lambda_small, R_lambda_small = [], []
    L_lambda_int, R_lambda_int = [], []
    L_lambda_large, R_lambda_large = [], []

    for index, row in df.iterrows():
        if index % 1000 == 0 and index > 0:
            click.echo("\t\tProcessed %d sites." % index)
        site_pos = row['site']

        flank_large_lower = site_pos - flank_large_len if site_pos - flank_large_len >= 0 else 0
        flank_large_upper = site_pos + flank_large_len if site_pos + flank_large_len < contig_length else contig_length - 1

        flank_int_lower = site_pos - flank_int_len if site_pos - flank_int_len >= 0 else 0
        flank_int_upper = site_pos + flank_int_len if site_pos + flank_int_len < contig_length else contig_length - 1

        flank_small_lower = site_pos - flank_small_len if site_pos - flank_small_len >= 0 else 0
        flank_small_upper = site_pos + flank_small_len if site_pos + flank_small_len < contig_length else contig_length - 1

        flanks_large = df.query("site <= {upper} & site >= {lower}".format(lower=flank_large_lower, upper=flank_large_upper))
        flanks_int = flanks_large.query("site <= {upper} & site >= {lower}".format(lower=flank_int_lower, upper=flank_int_upper))
        flanks_small = flanks_int.query("site <= {upper} & site >= {lower}".format(lower=flank_small_lower, upper=flank_small_upper))

        L_lambda_small.append(flanks_small.L_count.sum() / (flank_small_upper - flank_small_lower))
        R_lambda_small.append(flanks_small.R_count.sum() / (flank_small_upper - flank_small_lower))

        L_lambda_int.append(flanks_int.L_count.sum() / (flank_int_upper - flank_int_lower))
        R_lambda_int.append(flanks_int.R_count.sum() / (flank_int_upper - flank_int_lower))

        L_lambda_large.append(flanks_large.L_count.sum() / (flank_large_upper - flank_large_lower))
        R_lambda_large.append(flanks_large.R_count.sum() / (flank_large_upper - flank_large_lower))

    df.loc[:,'L_lambda_small'], df.loc[:,'R_lambda_small'] = L_lambda_small, R_lambda_small
    df.loc[:,'L_lambda_int'], df.loc[:,'R_lambda_int'] = L_lambda_int, R_lambda_int
    df.loc[:,'L_lambda_large'], df.loc[:,'R_lambda_large'] = L_lambda_large, R_lambda_large


    #print(df)

    # Apply the minimum count filter
    df.loc[:,'L_count'] = [c if c >= min_softclip_count else 0 for c in df['L_count']]
    df.loc[:,'R_count'] = [c if c >= min_softclip_count else 0 for c in df['R_count']]
    df = pd.DataFrame(df.query('L_count >= {min_softclip_count} | R_count >= {min_softclip_count}'.format(min_softclip_count=min_softclip_count)))

    # Choosing the maximum lambda

    df.loc[:,'L_lambda_max'] = df.apply(lambda r: max(r['L_lambda_bg'], r['L_lambda_small'], r['L_lambda_int'], r['L_lambda_large']), axis=1)
    df.loc[:,'R_lambda_max'] = df.apply(lambda r: max(r['R_lambda_bg'], r['R_lambda_small'], r['R_lambda_int'], r['R_lambda_large']), axis=1)

    df.loc[:,'L_pval'] = df.apply(lambda r: 1 - poisson.cdf(r['L_count'], r['L_lambda_max']) if r['L_count'] > 0 else np.nan, axis=1)
    df.loc[:,'R_pval'] = df.apply(lambda r: 1 - poisson.cdf(r['R_count'], r['R_lambda_max']) if r['R_count'] > 0 else np.nan, axis=1)


    df.loc[:,'L_pass'] = df.apply(lambda r: True if r['L_pval'] < pval_cutoff else False, axis=1)
    df.loc[:,'R_pass'] = df.apply(lambda r: True if r['R_pval'] < pval_cutoff else False, axis=1)

    # Now filter out everything that failed poisson cutoff.
    df = df.query("L_pass == True | R_pass == True")

    # Now remove sites that have no oppositely-oriented passing site nearby.
    df = misc.keep_sites_with_nearby_mates1(df, min_softclip_pair_distance, max_softclip_pair_distance)

    df = df.loc[:, ['contig', 'site', 'L_count', 'R_count', 'L_pval', 'R_pval', 'L_pass', 'R_pass']]
    df = df.sort_values(['contig', 'site'])
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
    return read.reference_name, read.get_reference_positions()[0] - 1

def right_softclip_site(read):
    return read.reference_name, read.get_reference_positions()[-1] + 1