import sys
from collections import defaultdict, OrderedDict
import pandas as pd
from mustache import alignment_tools

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

        out_df['contig'].append(contig)
        out_df['left_site'].append(left)
        out_df['right_site'].append(right)
        out_df['count'].append(count)
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





