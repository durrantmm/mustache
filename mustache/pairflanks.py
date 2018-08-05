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
from os.path import basename, join
from mustache.config import TMPDIR

verbose=True
logger = gogo.Gogo(__name__, verbose=verbose).logger

def get_flank_pairs(flanks, minimum_partner_runthrough_count, tmp_output_prefix=None, max_direct_repeat_length=21,
                    truncated_flank_length=40, ir_distance_from_end=15):

    logger.info("Finding all flank pairs within %d bases of each other ..." % max_direct_repeat_length)
    pairs = pair_all_nearby_flanks(flanks, max_direct_repeat_length)
    logger.info("Finding all inverted repeats in %d candidate pairs..." % pairs.shape[0])
    pairs = check_pairs_for_ir(pairs, truncated_flank_length, ir_distance_from_end, tmp_output_prefix)
    logger.info("Filtering pairs according to existence of inverted repeats, flank length difference, and read count difference...")
    filtered_pairs = filter_pairs(pairs, minimum_partner_runthrough_count)

    return filtered_pairs


def pair_all_nearby_flanks(flanks, max_direct_repeat_length):

    column_names = ['contig', 'index_5p', 'index_3p', 'pos_5p', 'pos_3p', 'softclip_count_5p', 'softclip_count_3p',
                    'runthrough_count_5p', 'runthrough_count_3p', 'seq_5p', 'seq_3p']
    outpairs = dict()

    for index, row in flanks.iterrows():

        if row.orient != 'R':
            continue

        contig, pos_3p = row['contig'], row['pos']
        softclip_count_3p, runthrough_count_3p, seq_3p = row['softclip_count'], row['runthrough_count'], row['consensus_seq']
        min_pos = pos_3p - max_direct_repeat_length

        candidate_pairs = flanks.query('contig == @contig & pos >= @min_pos & pos < @pos_3p & orient == "L"')

        for index2, row2 in candidate_pairs.iterrows():
            pos_5p, softclip_count_5p, runthrough_count_5p, seq_5p = row2['pos'], row2['softclip_count'], row2['runthrough_count'], row2['consensus_seq']
            outpairs[len(outpairs)] = [contig, index2, index, pos_5p, pos_3p, softclip_count_5p, softclip_count_3p,
                                       runthrough_count_3p, runthrough_count_5p, seq_5p, seq_3p]

    outpairs = pd.DataFrame.from_dict(outpairs, orient='index', columns=column_names)
    return outpairs


def check_pairs_for_ir(pairs, truncated_flank_length, ir_distance_from_end, tmp_output_prefix='mustache'):

    has_ir_all = []
    ir_5p_all = []
    ir_3p_all = []

    for index, row in pairs.iterrows():
        contig, index_5p, index_3p, pos_5p, pos_3p, softclip_count_5p, softclip_count_3p, \
        runthrough_count_5p, runthrough_count_3p, seq_5p, seq_3p = row

        trunc_seq_5p = truncate_sequence(seq_5p, truncated_flank_length, orient='R')
        trunc_seq_3p = truncate_sequence(seq_3p, truncated_flank_length, orient='L')

        combined_seq = str('N'*20).join([trunc_seq_5p, trunc_seq_3p])


        tmp_fasta_path = join(TMPDIR, tmp_output_prefix + '.' + str(randint(0, 1e20)) + '.fasta')
        fastatools.write_sequences_to_fasta([combined_seq], tmp_fasta_path)

        tmp_einverted_outfile = join(TMPDIR, tmp_output_prefix + '.' + str(randint(0, 1e20)) + '.out')
        tmp_einverted_outseq = join(TMPDIR, tmp_output_prefix + '.' + str(randint(0, 1e20)) + '.fa')
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

def filter_pairs(pairs, minimum_partner_runthrough_count):
    pairs.loc[:, 'direct_repeat_length'] = pairs.loc[:, 'pos_3p'] - pairs.loc[:, 'pos_5p'] -1
    pairs.query("direct_repeat_length == 0 | (runthrough_count_5p >= @minimum_partner_runthrough_count & "
                "runthrough_count_3p >= @minimum_partner_runthrough_count)", inplace=True)

    pairs['IR_length'] = np.array([len(seq) if seq is not None else 0 for seq in list(pairs['IR_5p'])])
    pairs['difflength'] = abs(np.array(list(map(len, pairs['seq_5p']))) - np.array(list(map(len, pairs['seq_3p']))))
    pairs['diffcount'] = abs(pairs['softclip_count_5p'] - pairs['softclip_count_3p'])
    pairs['ignore_pair'] = False
    sorted_pairs = pairs.sort_values(['IR_length', 'difflength', 'diffcount'], ascending=[False, True, True])

    for index, row in sorted_pairs.iterrows():
        if not row['ignore_pair']:
            index_5p, index_3p = row['index_5p'], row['index_3p']
            has_pair_member = sorted_pairs.query("index_5p == @index_5p | index_3p == @index_3p")
            has_pair_member_indices = has_pair_member.query("index_5p != @index_5p | index_3p != @index_3p").index.values
            sorted_pairs.loc[has_pair_member_indices, 'ignore_pair'] = True

    filtered_pairs = sorted_pairs.query("ignore_pair == False").loc[:,
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

def _pairflanks(flanksfile, minimum_partner_runthrough_count, output_file):
    tmp_output_prefix = '.'.join(basename(output_file).split('.')[:-1])

    flanks = pd.read_csv(flanksfile, sep='\t')

    if flanks.shape[0] == 0:
        logger.info("No flanks found in the input file...")
        if output_file:
            flanks.to_csv(output_file, sep='\t', index=False)
        return flanks

    else:

        flank_pairs = get_flank_pairs(flanks, minimum_partner_runthrough_count, tmp_output_prefix=tmp_output_prefix)

        logger.info("Identified %d flank pairs with inverted repeats..." % flank_pairs.shape[0])
        if output_file:
            logger.info("Saving results to file %s" % output_file)
            flank_pairs.to_csv(output_file, sep='\t', index=False)

        return flank_pairs

@click.command()
@click.argument('flanksfile', type=click.Path(exists=True))
@click.option('--output_file', '-o', default=None, help="The output file to save the results.")
@click.option('--minimum_partner_runthrough_count', '-minrtcount', default=5,
              help="If two flanks are more than one base apart, then you would expect each softclipped site to"
                   "also show signs of non-softclippd reads, signifying a direct repeat. This is the number of"
                   "such runthrough reads required for a pair to be considered valid.")
def pairflanks(flanksfile, minimum_partner_runthrough_count, output_file=None):
    if not output_file:
        output_file = '.'.join(['mustache', basename(flanksfile).split('.')[0], 'pairflanks.tsv'])
    _pairflanks(flanksfile, minimum_partner_runthrough_count, output_file)


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