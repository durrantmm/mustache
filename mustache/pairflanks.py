import warnings
warnings.filterwarnings("ignore")
import sys
import click
import pygogo as gogo
import pandas as pd
import numpy as np
from snakemake import shell
from random import randint
from mustache import fastatools, embosstools
from os.path import basename, join, dirname
from Bio import SeqIO
import pysam
from collections import defaultdict
verbose=True
logger = gogo.Gogo(__name__, verbose=verbose).logger

def get_flank_pairs(flanks, tmp_dir, tmp_output_prefix=None, max_direct_repeat_length=20,
                    truncated_flank_length=40, ir_distance_from_end=15):

    logger.info("Finding all flank pairs within %d bases of each other ..." % max_direct_repeat_length)
    pairs = pair_all_nearby_flanks(flanks, max_direct_repeat_length)
    logger.info("Finding all inverted repeats at termini in %d candidate pairs..." % pairs.shape[0])
    pairs = check_pairs_for_ir(pairs, truncated_flank_length, ir_distance_from_end, tmp_dir, tmp_output_prefix)
    logger.info("Filtering pairs according to existence of inverted repeats, flank length difference, and read count difference...")
    filtered_pairs = filter_pairs(pairs)

    return filtered_pairs


def pair_all_nearby_flanks(flanks, max_direct_repeat_length):

    column_names = ['contig', 'index_5p', 'index_3p', 'pos_5p', 'pos_3p', 'softclip_count_5p', 'softclip_count_3p',
                    'runthrough_count_5p', 'runthrough_count_3p']

    column_names += ['seq_5p', 'seq_3p']

    outpairs = dict()

    for index_5p, row in flanks.iterrows():

        if row.orient != 'R':
            continue

        contig, pos_5p = row['contig'], row['pos']
        softclip_count_5p, runthrough_count_5p, seq_5p = row['softclip_count'], row['runthrough_count'], row['consensus_seq']

        min_pos = pos_5p - max_direct_repeat_length

        candidate_pairs = flanks.query('contig == @contig & pos > @min_pos & pos < @pos_5p & orient == "L"')

        for index_3p, row2 in candidate_pairs.iterrows():
            pos_3p, softclip_count_3p, runthrough_count_3p, seq_3p = row2['pos'], row2['softclip_count'], \
                                                                     row2['runthrough_count'], row2['consensus_seq']

            out = [contig, index_5p, index_3p, pos_5p, pos_3p, softclip_count_5p, softclip_count_3p, runthrough_count_5p, runthrough_count_3p]



            out += [seq_5p, seq_3p]
            outpairs[len(outpairs)] =  out

    outpairs = pd.DataFrame.from_dict(outpairs, orient='index', columns=column_names)
    return outpairs


def check_pairs_for_ir(pairs, truncated_flank_length, ir_distance_from_end, tmp_dir, tmp_output_prefix='mustache'):

    has_ir_all = []
    ir_5p_all = []
    ir_3p_all = []

    for index, row in pairs.iterrows():

        contig, index_5p, index_3p, pos_5p, pos_3p, softclip_count_5p, softclip_count_3p, \
        runthrough_count_5p, runthrough_count_3p, seq_5p, seq_3p = row

        trunc_seq_5p = truncate_sequence(seq_5p, truncated_flank_length, orient='R')
        trunc_seq_3p = truncate_sequence(seq_3p, truncated_flank_length, orient='L')

        combined_seq = str('N'*20).join([trunc_seq_5p, trunc_seq_3p])


        tmp_fasta_path = join(tmp_dir, tmp_output_prefix + '.' + str(randint(0, 1e20)) + '.fasta')
        fastatools.write_sequences_to_fasta([combined_seq], tmp_fasta_path)

        tmp_einverted_outfile = join(tmp_dir, tmp_output_prefix + '.' + str(randint(0, 1e20)) + '.out')
        tmp_einverted_outseq = join(tmp_dir, tmp_output_prefix + '.' + str(randint(0, 1e20)) + '.fa')
        embosstools.run_einverted(tmp_fasta_path, outfile=tmp_einverted_outfile, outseq=tmp_einverted_outseq)

        has_ir = False
        ir_length = 0
        keep_ir1 = None
        keep_ir2 = None
        for ir1, ir2 in embosstools.read_emboss_seq_results(tmp_einverted_outseq):
            if pair_has_ir(ir1, ir2, ir_distance_from_end, len(combined_seq)):
                has_ir = True
                if len(ir1.seq) > ir_length:
                    keep_ir1 = ir1.seq
                    keep_ir2 = ir2.seq

        shell('rm -f %s' % tmp_fasta_path)
        shell('rm -f %s' % tmp_einverted_outfile)
        shell('rm -f %s' % tmp_einverted_outseq)

        has_ir_all.append(has_ir)
        ir_3p_all.append(keep_ir1)
        ir_5p_all.append(keep_ir2)

    pairs['has_IR'] = has_ir_all
    pairs['IR_5p'] = ir_5p_all
    pairs['IR_3p'] = ir_3p_all

    return pairs


def filter_pairs(pairs):
    pairs.loc[:, 'direct_repeat_length'] = pairs.loc[:, 'pos_5p'] - pairs.loc[:, 'pos_3p'] -1

    pairs['IR_length'] = np.array([len(seq) if seq is not None else 0 for seq in list(pairs.loc[:,'IR_5p'])])
    pairs['difflength'] = abs(np.array(list(map(len, pairs['seq_5p']))) - np.array(list(map(len, pairs['seq_3p']))))
    pairs['diffcount'] = abs(pairs['softclip_count_5p'] - pairs['softclip_count_3p'])
    pairs['ignore_pair'] = False
    sorted_pairs = pairs.sort_values(['IR_length', 'difflength', 'diffcount'], ascending=[False, True, True])

    keep_row = []

    for index in sorted_pairs.index:
        row = sorted_pairs.loc[index, :]
        if row['pos_3p'] == 15377:
            print(row)

        if not row['ignore_pair']:
            keep_row.append(True)
            index_5p, index_3p = row['index_5p'], row['index_3p']
            has_pair_member = sorted_pairs.query("index_5p == @index_5p | index_3p == @index_3p")
            has_pair_member_indices = has_pair_member.query("index_5p != @index_5p | index_3p != @index_3p").index.values
            sorted_pairs.loc[has_pair_member_indices, 'ignore_pair'] = True
        else:
            keep_row.append(False)

    sorted_pairs.loc[:, 'keep_pair'] = np.array(keep_row)
    filtered_pairs = sorted_pairs.query("keep_pair == True").loc[:,
                     ['contig', 'pos_5p', 'pos_3p', 'softclip_count_5p', 'softclip_count_3p', 'runthrough_count_5p',
                      'runthrough_count_3p', 'has_IR', 'IR_length', 'IR_5p', 'IR_3p',
                      'seq_5p', 'seq_3p']].sort_values(['contig', 'pos_5p', 'pos_3p'])

    return filtered_pairs

def pair_has_ir(ir1, ir2, ir_distance_from_end, seqlen):
    if ir_near_5prime_end(ir1, ir_distance_from_end) and ir_near_3prime_end(ir2, ir_distance_from_end, seqlen):
        return True
    return False

def ir_near_5prime_end(ir1, ir_distance_from_end):
    if ir1.ir_pos_5p <= ir_distance_from_end:
        return True
    return False

def ir_near_3prime_end(ir2, ir_distance_from_end, seqlen):
    if ir2.ir_pos_3p >= (seqlen - ir_distance_from_end):
        return True
    return False

def truncate_sequence(seq, truncated_seq_length, orient='R'):
    truncated_seq = seq
    if len(truncated_seq) > truncated_seq_length:
        if orient == 'R':
            truncated_seq = truncated_seq[:truncated_seq_length]
        elif orient == 'L':
            truncated_seq = truncated_seq[-truncated_seq_length:]
    return truncated_seq

def get_direct_repeats(flank_pairs, bamfile, genome):

    genome_dict = {rec.id: rec.seq for rec in SeqIO.parse(genome, 'fasta')}
    positions = get_reference_direct_repeats(flank_pairs, genome_dict)
    positions = get_read_direct_repeats(positions, genome_dict, bamfile)

    flank_pairs = flank_pairs.merge(positions, how='left')
    return flank_pairs

def get_read_direct_repeats(positions, genome_dict, bamfile, target_region_size=50):

    bam = pysam.AlignmentFile(bamfile, 'rb')

    direct_repeats = []
    for index, row in positions.iterrows():
        contig, start, end = row['contig'], row['pos_3p'], row['pos_5p']

        direct_repeat_center = round((end + start) / 2)
        expanded_start = int(direct_repeat_center - (target_region_size / 2))
        expanded_end = int(direct_repeat_center + (target_region_size / 2))

        add_start_n = 0
        add_end_n = 0
        if expanded_start < 0:
            add_start_n = abs(expanded_start)
            expanded_start = 0
        if expanded_end > len(genome_dict[contig]):
            add_end_n = expanded_end - len(genome_dict[contig])
            expanded_end = len(genome_dict[contig])

        target_region = 'N' * add_start_n + genome_dict[contig][expanded_start:expanded_end] + 'N' * add_end_n
        target_region_positions = range(expanded_start, expanded_end)

        target_region_reads = initialize_target_region_reads(target_region, expanded_start, expanded_end)
        for read in bam.fetch(contig, expanded_start, expanded_end):
            ref_positions = read.get_reference_positions(full_length=True)
            read_query = read.query_sequence
            read_qualities = read.query_qualities

            i = 0
            for pos in ref_positions:
                if pos is not None and pos in target_region_reads:
                    target_region_reads[pos][read_query[i]] += read_qualities[i]
                i += 1

        consensus_target_region = get_consensus_target_region(target_region_reads)
        consensus_direct_repeat = consensus_target_region[target_region_positions.index((start+1)):target_region_positions.index(end)]

        direct_repeats.append(consensus_direct_repeat)

    positions['direct_repeat_reads_consensus'] = direct_repeats

    return positions

def get_consensus_target_region(target_region_reads):
    start, end = min(target_region_reads.keys()), max(target_region_reads.keys())+1
    consensus = ''
    for pos in range(start, end):
        best_qual = 0
        best_nuc = ''
        for nuc in target_region_reads[pos]:
            qual = target_region_reads[pos][nuc]
            if qual > best_qual:
                best_qual = qual
                best_nuc = nuc
        consensus += best_nuc
    return consensus


def initialize_target_region_reads(target_region, expanded_start, expanded_end):
    target_region_reads = defaultdict(lambda: defaultdict(int))
    i = 0
    for pos in range(expanded_start, expanded_end):
        target_region_reads[pos][target_region[i]] += 1
        i += 1
    return target_region_reads


def get_reference_direct_repeats(flank_pairs, genome_dict, target_region_size=50):
    positions = flank_pairs.loc[:, ['contig', 'pos_5p', 'pos_3p']].drop_duplicates().reset_index(drop=True)

    direct_repeats = []
    for index, row in positions.iterrows():
        contig, start, end = row['contig'], row['pos_3p'], row['pos_5p']

        direct_repeat = genome_dict[contig][(start+1):end]
        direct_repeats.append(''.join(direct_repeat))


    positions['direct_repeat_reference'] = direct_repeats
    return positions


def _pairflanks(flanksfile, bamfile, genome, max_direct_repeat_length, output_file):
    tmp_output_prefix = '.'.join(basename(output_file).split('.')[:-1])

    flanks = pd.read_csv(flanksfile, sep='\t')

    if flanks.shape[0] == 0:
        logger.info("No flanks found in the input file...")
        flank_pairs = pd.DataFrame(
            columns=['pair_id', 'contig', 'pos_5p', 'pos_3p', 'softclip_count_5p', 'softclip_count_3p',
                     'runthrough_count_5p', 'runthrough_count_3p', 'has_IR', 'IR_length',
                     'IR_5p', 'IR_3p', 'seq_5p', 'seq_3p', 'direct_repeat_reference',
                     'direct_repeat_reads_consensus'])

        if output_file:
            flank_pairs.to_csv(output_file, sep='\t', index=False)
        return flank_pairs

    else:

        flank_pairs = get_flank_pairs(flanks, dirname(output_file), tmp_output_prefix=tmp_output_prefix,
                                      max_direct_repeat_length=max_direct_repeat_length)
        logger.info("Identified %d flank pairs in total..." % flank_pairs.shape[0])
        logger.info("Identified %d flank pairs with inverted repeats..." % flank_pairs.query('has_IR==True').shape[0])
        logger.info("Getting direct repeats and surrounding genomic region...")
        flank_pairs = get_direct_repeats(flank_pairs, bamfile, genome)

        if output_file:
            logger.info("Saving results to file %s" % output_file)
            flank_pairs.reset_index(drop=True)
            flank_pairs.index = flank_pairs.index + 1
            flank_pairs.index.name = 'pair_id'
            flank_pairs.to_csv(output_file, sep='\t')

        return flank_pairs