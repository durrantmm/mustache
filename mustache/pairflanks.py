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

verbose=True
logger = gogo.Gogo(__name__, verbose=verbose).logger

def get_flank_pairs(flanks, tmp_dir, tmp_output_prefix=None, max_direct_repeat_length=18,
                    truncated_flank_length=40, ir_distance_from_end=15):

    logger.info("Finding all flank pairs within %d bases of each other ..." % max_direct_repeat_length)
    pairs = pair_all_nearby_flanks(flanks, max_direct_repeat_length)
    logger.info("Finding all inverted repeats in %d candidate pairs..." % pairs.shape[0])
    pairs = check_pairs_for_ir(pairs, truncated_flank_length, ir_distance_from_end, tmp_dir, tmp_output_prefix)
    logger.info("Filtering pairs according to existence of inverted repeats, flank length difference, and read count difference...")
    filtered_pairs = filter_pairs(pairs)

    return filtered_pairs


def pair_all_nearby_flanks(flanks, max_direct_repeat_length):

    column_names = ['contig', 'index_5p', 'index_3p', 'pos_5p', 'pos_3p', 'softclip_count_5p', 'softclip_count_3p',
                    'runthrough_count_5p', 'runthrough_count_3p']

    if 'extended' in list(flanks.loc[0,:].keys()):
        column_names += ['extended_5p', 'extended_3p']

    column_names += ['seq_5p', 'seq_3p']

    outpairs = dict()

    for index_5p, row in flanks.iterrows():

        if row.orient != 'R':
            continue

        contig, pos_5p = row['contig'], row['pos']
        softclip_count_5p, runthrough_count_5p, seq_5p = row['softclip_count'], row['runthrough_count'], row['consensus_seq']

        if 'extended' in list(row.keys()):
            extended_5p = row['extended']

        min_pos = pos_5p - max_direct_repeat_length

        candidate_pairs = flanks.query('contig == @contig & pos > @min_pos & pos < @pos_5p & orient == "L"')

        for index_3p, row2 in candidate_pairs.iterrows():
            pos_3p, softclip_count_3p, runthrough_count_3p, seq_3p = row2['pos'], row2['softclip_count'], \
                                                                     row2['runthrough_count'], row2['consensus_seq']
            if 'extended' in list(row2.keys()):
                extended_3p = row2['extended']

            out = [contig, index_5p, index_3p, pos_5p, pos_3p, softclip_count_5p, softclip_count_3p, runthrough_count_5p, runthrough_count_3p]

            if 'extended' in list(row2.keys()):
                out += [extended_5p, extended_3p]

            out += [seq_5p, seq_3p]
            outpairs[len(outpairs)] =  out

    outpairs = pd.DataFrame.from_dict(outpairs, orient='index', columns=column_names)
    return outpairs


def check_pairs_for_ir(pairs, truncated_flank_length, ir_distance_from_end, tmp_dir, tmp_output_prefix='mustache', ):

    has_ir_all = []
    ir_5p_all = []
    ir_3p_all = []

    for index, row in pairs.iterrows():
        if 'extended_5p' in list(row.keys()):
            contig, index_5p, index_3p, pos_5p, pos_3p, softclip_count_5p, softclip_count_3p, \
            runthrough_count_5p, runthrough_count_3p, extended_5p, extended_3p, seq_5p, seq_3p = row
        else:
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
                      'runthrough_count_3p', 'extended_5p', 'extended_3p', 'has_IR', 'IR_length', 'IR_5p', 'IR_3p',
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

def _pairflanks(flanksfile, output_file):
    tmp_output_prefix = '.'.join(basename(output_file).split('.')[:-1])

    flanks = pd.read_csv(flanksfile, sep='\t')

    if flanks.shape[0] == 0:
        logger.info("No flanks found in the input file...")
        if output_file:
            flanks.to_csv(output_file, sep='\t', index=False)
        return flanks

    else:

        flank_pairs = get_flank_pairs(flanks, dirname(output_file), tmp_output_prefix=tmp_output_prefix)

        logger.info("Identified %d flank pairs with inverted repeats..." % flank_pairs.shape[0])
        if output_file:
            logger.info("Saving results to file %s" % output_file)
            flank_pairs.to_csv(output_file, sep='\t', index=False)

        return flank_pairs

@click.command()
@click.argument('flanksfile', type=click.Path(exists=True))
@click.option('--output_file', '-o', default='mustache.pairflanks.tsv', help="The output file to save the results.")
def pairflanks(flanksfile, output_file=None):
    _pairflanks(flanksfile, output_file)


if __name__ == '__main__':
    pairflanks()
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