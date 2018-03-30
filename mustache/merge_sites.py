import sys
from collections import defaultdict, OrderedDict
import pandas as pd
from mustache import alignment_tools
import numpy as np

def merge(sites_files):

    sites = defaultdict(lambda: defaultdict(list))

    for f in sites_files:
        df = pd.read_csv(f, sep='\t')
        df = df[['contig', 'left_site', 'right_site', 'right_assembly', 'left_assembly', 'merged_assembly']]
        for row in range(len(df)):

            contig = df.iloc[row,]['contig']
            left_site = df.iloc[row,]['left_site']
            right_site = df.iloc[row,]['right_site']

            right_assembly = df.iloc[row,]['right_assembly']
            left_assembly = df.iloc[row,]['left_assembly']
            merged_assembly = df.iloc[row,]['merged_assembly']

            sites[(contig, left_site, right_site)]['right'].append(right_assembly)
            sites[(contig, left_site, right_site)]['left'].append(left_assembly)
            sites[(contig, left_site, right_site)]['merged'].append(merged_assembly)


    out_df = OrderedDict([('contig', []), ('left_site',[]), ('right_site',[]), ('count',[]),
                          ('right_assembly',[]), ('left_assembly',[]), ('merged_assembly',[])])
    for s in sites:

        contig, left, right = s
        right_assemblies, left_assemblies, merged_assemblies = sites[s]['right'], sites[s]['left'], sites[s]['merged']

        right_asm, left_asm, merged_asm, count = merge_assemblies(right_assemblies, left_assemblies, merged_assemblies)

        right_length, left_length = len(right_asm), len(left_asm)
        right_repeats, left_repeats = has_repeats(right_asm), has_repeats(left_asm)


        merged_length, merged_repeats = np.nan, np.nan

        if not np.isnan(merged_asm):
            merged_length, merged_repeats = len(merged_asm), has_repeats(merged_asm)

        out_df['contig'].append(contig)
        out_df['left_site'].append(left)
        out_df['right_site'].append(right)
        out_df['count'].append(count)
        out_df['right_length'].append(right_length)
        out_df['left_length'].append(left_length)
        out_df['merged_length'].append(merged_length)
        out_df['right_repeats'].append(right_repeats)
        out_df['left_repeats'].append(left_repeats)
        out_df['merged_repeats'].append(merged_length)
        out_df['right_assembly'].append(right_asm)
        out_df['left_assembly'].append(left_asm)
        out_df['merged_assembly'].append(merged_asm)

    out_df = pd.DataFrame(out_df)

    return out_df

def merge_assemblies(right_assemblies, left_assemblies, merged_assemblies):

    if len(left_assemblies) != len(right_assemblies):
        print("Error: Must have same number of left and right assemblies.")
        sys.exit()

    if len(left_assemblies) == 1:
        return right_assemblies[0], left_assemblies[0], merged_assemblies[0], 1

    merged_right_assembly = right_assemblies[0]
    for asm in right_assemblies[1:]:
        length = min([len(merged_right_assembly), len(asm)])
        score = alignment_tools.left_alignment_score(merged_right_assembly, asm)

        if score == 0:
            continue
        elif length / float(score) > 0.9:
            merged_right_assembly = max([merged_right_assembly, asm], key=len)

    merged_left_assembly = left_assemblies[0]
    for asm in left_assemblies[1:]:
        length = min([len(merged_left_assembly), len(asm)])
        score = alignment_tools.right_alignment_score(merged_left_assembly, asm)

        if score == 0:
            continue
        elif length / float(score) > 0.9:
            merged_left_assembly = max([merged_left_assembly, asm], key=len)

    l, r, final_merged_assembly = alignment_tools.merge_flank_assemblies(merged_right_assembly, merged_left_assembly, verbose=False)

    return merged_right_assembly, merged_left_assembly, final_merged_assembly, len(right_assemblies)

def second_largest(numbers):
    num_max = max(numbers)
    found_max = False
    second_largest = 0
    for n in numbers:
        if n == num_max and not found_max:
            found_max = True
        elif n == num_max and found_max:
            second_largest = n
        elif n > second_largest:
            second_largest = n
    return second_largest

def has_repeats(dna, ratio_cutoff=0.1):
    kmer_lengths = list()
    ratios = list()

    for k in range(1, len(dna)):

        kmer_ranges = range(0, len(dna), k)

        if len(list(kmer_ranges)) == 2:
            break

        kmer_lengths.append(k)
        kmer_dict = defaultdict(int)
        for start in range(0, len(dna), k):
            if start + k > len(dna):
                break
            kmer = dna[start:(start + k)]
            kmer_dict[kmer] += 1

        max_count = max(kmer_dict.values())
        second_max = second_largest(kmer_dict.values())

        ratio = second_max / max_count
        ratios.append(ratio)

    if min(ratios) < ratio_cutoff:
        return True
    else:
        return False



