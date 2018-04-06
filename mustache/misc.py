import sys
import numpy as np
from scipy.stats import norm
import pandas as pd

def calculate_flank_length(bam_file, max_mapping_distance, flank_length_percentile, max_reads = 100000):
    fragment_lengths = []
    count = 0
    for read in bam_file.fetch():
        if count > max_reads:
            break
        if not read.is_reverse:
            if read.tlen != 0 and abs(read.tlen) < max_mapping_distance:
                fragment_length = abs(read.tlen)
                fragment_lengths.append(fragment_length)
                count += 1
    mean = np.mean(fragment_lengths)
    std = np.std(fragment_lengths)

    if std == 0:
        flank_length = mean
    else:
        flank_length = int(norm.ppf(flank_length_percentile, mean, std))

    return flank_length


def calculate_average_read_length(bam_file, max_reads = 100000):
    read_lengths = []
    count = 0
    for read in bam_file.fetch():
        count += 1
        if count > max_reads:
            break
        read_lengths.append(len(read.query_sequence))

    average = int(np.mean(read_lengths))

    return average


def get_unique_kmers(seqs, k):
    kmers = set()

    for seq in seqs:
        if len(seq) < k:
            continue
        for i in range(k, len(seq)):
            kmers.add(seq[(i - k):i])

    return sorted(list(kmers))


def calculate_maximum_softclip_length(bam_file, max_reads = 100000):
    max_softclip_length = 0
    count = 0
    for read in bam_file.fetch():
        if count > max_reads:
            break

        softclip_length = 0
        if read.cigartuples[0][0] == 4 and read.cigartuples[-1][0] != 4:
            softclip_length = read.cigartuples[0][1]
        elif read.cigartuples[0][0] != 4 and read.cigartuples[-1][0] == 4:
            softclip_length = read.cigartuples[-1][1]

        if softclip_length != 0:
            count += 1

        if softclip_length > max_softclip_length:
            max_softclip_length = softclip_length


    return max_softclip_length


def calculate_maximum_read_length(bam_file, max_reads=100000):
    max_read_length = 0
    count = 0
    for read in bam_file.fetch():

        count += 1

        if count > max_reads:
            break

        if len(read.query_sequence) > max_read_length:
            max_read_length = len(read.query_sequence)

    return max_read_length


def revcomp(read):
    reversed_seq = ''
    for l in reversed(read):
        if l == 'A':
            reversed_seq += 'T'
        elif l == 'T':
            reversed_seq += 'A'
        elif l == 'C':
            reversed_seq += 'G'
        elif l == 'G':
            reversed_seq += 'C'
        else:
            reversed_seq += l
    return reversed_seq


def wrap_string(string, wrap_cutoff=70, newline_char=''):
    out = ''
    count = 0
    for c in string:
        if len(out) == 0:
            out += newline_char
            out += c
        elif count == wrap_cutoff:
            out += '\n'
            out += newline_char
            out += c
            count = 0
        else:
            out += c
        count += 1

    return out


def get_bam_contig_dict(bam_file):
    contig_dict = {}
    for contig in bam_file.header['SQ']:
        contig_dict[contig['SN']] = contig['LN']
    return contig_dict


def keep_sites_with_nearby_mates1(df, min_distance, max_distance, contig_column='contig', position_column='site'):
    L_has_nearby_mate = []
    R_has_nearby_mate = []
    for index, row in df.iterrows():
        site_pos = row[position_column]
        site_contig = row[contig_column]
        L_has_mate, R_has_mate = False, False
        if row['L_pass']:
            nearby_mates = df.query(
                'R_pass == True & {position_column} - {site_pos} >= {min_softclip_pair_distance} & {position_column} - {site_pos} <= {max_softclip_pair_distance} & contig == "{site_contig}"'.format(
                    site_pos=site_pos, min_softclip_pair_distance=min_distance,
                    max_softclip_pair_distance=max_distance,
                    position_column=position_column, site_contig=site_contig))
            if nearby_mates.shape[0] > 0:
                L_has_mate = True

        if row['R_pass']:
            nearby_mates = df.query(
                'L_pass == True & {site_pos} - {position_column} >= {min_softclip_pair_distance} & {site_pos} - {position_column} <= {max_softclip_pair_distance} & contig == "{site_contig}"'.format(
                    site_pos=site_pos, min_softclip_pair_distance=min_distance,
                    max_softclip_pair_distance=max_distance,
                    position_column=position_column, site_contig=site_contig))
            if nearby_mates.shape[0] > 0:
                R_has_mate = True

        L_has_nearby_mate.append(L_has_mate)
        R_has_nearby_mate.append(R_has_mate)

    df.loc[:, 'L_has_nearby_mate'] = L_has_nearby_mate
    df.loc[:, 'R_has_nearby_mate'] = R_has_nearby_mate

    out_df = pd.DataFrame(df.query("L_has_nearby_mate == True | R_has_nearby_mate == True"))
    out_df.drop(['L_has_nearby_mate', 'R_has_nearby_mate'], axis=1, inplace=True)
    return out_df

def keep_sites_with_nearby_mates2(df, min_distance, max_distance, contig_column='contig', position_column='site'):
    L_has_nearby_mate = []
    R_has_nearby_mate = []
    for index, row in df.iterrows():
        site_contig = row[contig_column]
        site_pos = row[position_column]

        L_has_mate, R_has_mate = False, False
        if row['orientation'] == "L":
            nearby_mates = df.query(
                'orientation == "R" & {position_column} - {site_pos} >= {min_softclip_pair_distance} & {position_column} - {site_pos} <= {max_softclip_pair_distance} & contig == "{site_contig}"'.format(
                    site_pos=site_pos, min_softclip_pair_distance=min_distance,
                    max_softclip_pair_distance=max_distance,
                    position_column=position_column, site_contig=site_contig))
            if nearby_mates.shape[0] > 0:
                L_has_mate = True

        if row['orientation'] == "R":
            nearby_mates = df.query(
                'orientation == "L" & {site_pos} - {position_column} >= {min_softclip_pair_distance} & {site_pos} - {position_column} <= {max_softclip_pair_distance} & contig == "{site_contig}"'.format(
                    site_pos=site_pos, min_softclip_pair_distance=min_distance,
                    max_softclip_pair_distance=max_distance,
                    position_column=position_column, site_contig=site_contig))
            if nearby_mates.shape[0] > 0:
                R_has_mate = True

        L_has_nearby_mate.append(L_has_mate)
        R_has_nearby_mate.append(R_has_mate)

    df.loc[:, 'L_has_nearby_mate'] = L_has_nearby_mate
    df.loc[:, 'R_has_nearby_mate'] = R_has_nearby_mate

    out_df = pd.DataFrame(df.query("L_has_nearby_mate == True | R_has_nearby_mate == True"))
    out_df.drop(['L_has_nearby_mate', 'R_has_nearby_mate'], axis=1, inplace=True)
    return out_df


def get_nearby_sites(df, row, min_distance, max_distance, contig_column='contig', position_column='site'):

    site_pos = row[position_column]
    site_contig = row[contig_column]

    if row['orientation'] == 'L':
        nearby_mates = df.query(
            'orientation == "R" & {position_column} - {site_pos} >= {min_softclip_pair_distance} & {position_column} - {site_pos} <= {max_softclip_pair_distance} & contig == "{site_contig}"'.format(
                site_pos=site_pos, min_softclip_pair_distance=min_distance,
                max_softclip_pair_distance=max_distance,
                position_column=position_column, site_contig=site_contig))
    else:
        nearby_mates = df.query(
            'orientation == "L" & {site_pos} - {position_column} >= {min_softclip_pair_distance} & {site_pos} - {position_column} <= {max_softclip_pair_distance} & contig == "{site_contig}"'.format(
                site_pos=site_pos, min_softclip_pair_distance=min_distance,
                max_softclip_pair_distance=max_distance,
                position_column=position_column, site_contig=site_contig))

    return nearby_mates