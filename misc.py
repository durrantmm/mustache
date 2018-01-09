import numpy as np
import sys
from scipy.stats import norm

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