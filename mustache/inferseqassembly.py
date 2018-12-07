import sys
import warnings
warnings.filterwarnings("ignore")
import click
import pandas as pd
from mustache import fastatools
from mustache import bowtie2tools
from mustache import sctools
from mustache import misc
from mustache import pysamtools
import pygogo as gogo
import pysam
from Bio import SeqIO
from collections import OrderedDict
from os.path import dirname, join
from random import randint
from snakemake import shell
from collections import defaultdict

verbose=True
logger = gogo.Gogo(__name__, verbose=verbose).logger


import warnings
warnings.filterwarnings("ignore", message="numpy.dtype size changed")
warnings.filterwarnings("ignore", message="numpy.ufunc size changed")


def _inferseq_assembly(pairsfile, bamfile, inferseq_assembly, inferseq_reference, min_perc_identity,
                       max_internal_softclip_prop, max_inferseq_size, min_inferseq_size, keep_intermediate, output_file):

    index_genome(inferseq_assembly)

    bam = pysam.AlignmentFile(bamfile, 'rb')
    reference_genome_dict = {rec.id: rec.seq for rec in SeqIO.parse(inferseq_reference, 'fasta')}
    assembly_genome_dict = {rec.id: rec.seq for rec in SeqIO.parse(inferseq_assembly, 'fasta')}

    tmp_dir = dirname(output_file)

    pairs = pd.read_csv(pairsfile, sep='\t')
    handle_empty_pairsfile(pairs, output_file)

    logger.info("Aligning pairs to assembly...")
    assembly_flanks_fasta_prefix = write_flanks_to_align_to_assembly(pairs, bam, reference_genome_dict, tmp_dir)
    assembly_outbam = join(tmp_dir, 'mustache.inferseq_assembly.' + str(randint(0, 1e20)) + '.bam')
    bowtie2tools.align_fasta_to_genome(
        assembly_flanks_fasta_prefix+'.fasta',
        inferseq_assembly, assembly_outbam, silence=True,
        additional_flags='--all'
    )

    logger.info("Inferring sequences from pairs aligned to assembly and in context...")
    sequences_inferred_from_assembly_with_context = infer_sequences_with_context(assembly_outbam, assembly_genome_dict,
                                                                                 min_perc_identity, max_internal_softclip_prop)
    logger.info("Inferring sequences from pairs aligned to assembly...")
    sequences_inferred_from_assembly_without_context = infer_sequences_without_context(assembly_outbam, assembly_genome_dict,
                                                                                       min_perc_identity, max_internal_softclip_prop)

    if not keep_intermediate:
        shell('rm {fasta_prefix}* {outbam}*'.format(fasta_prefix=assembly_flanks_fasta_prefix, outbam=assembly_outbam))

    method1 = make_dataframe(sequences_inferred_from_assembly_with_context, method='inferred_assembly_with_context')
    method2 = make_dataframe(sequences_inferred_from_assembly_without_context, method='inferred_assembly_without_context')

    all_inferred_results = method1.append(method2, ignore_index=True).sort_values(by=['pair_id', 'method'])

    all_inferred_results.loc[:, 'pair_id'] = list(map(str, map(int, list(all_inferred_results['pair_id']))))
    all_inferred_results = all_inferred_results.query("inferred_seq_length >= @min_inferseq_size")
    all_inferred_results = all_inferred_results.query("inferred_seq_length <= @max_inferseq_size")

    logger.info("Writing results to file %s..." % output_file)
    all_inferred_results.to_csv(output_file, sep='\t', index=False)


def make_dataframe(inferred_sequences, method=None):
    outdict = OrderedDict([("pair_id", []), ("method", []), ("loc", []),
                           ("inferred_seq_length", []), ("inferred_seq", [])])
    for pair_id in inferred_sequences:
        for result in inferred_sequences[pair_id]:
            outdict['pair_id'].append(int(pair_id.split('_')[0]))
            outdict['method'].append(method)
            outdict['loc'].append(str(result[0]))
            outdict['inferred_seq_length'].append(result[1])
            outdict['inferred_seq'].append(''.join(result[2]))

    outdf = pd.DataFrame.from_dict(outdict)

    return outdf


def infer_sequences_with_context(bam_file, genome_dict, min_perc_identity, max_internal_softclip_prop):
    bam = pysam.AlignmentFile(bam_file, 'rb')
    keep_reads = prefilter_context_reads(bam, genome_dict, min_perc_identity)
    keep_pairs = get_pairs(keep_reads)
    keep_pairs = filter_context_pairs(keep_pairs, max_internal_softclip_prop)

    for pair_id in keep_pairs:
        keep_pairs[pair_id] = keep_best_alignment_score(keep_pairs[pair_id])

    inferred_sequences = defaultdict(list)
    for pair_id in keep_pairs:
        inferred_sequences[pair_id] = get_inferred_sequences(keep_pairs[pair_id], genome_dict, add_softclipped_bases=True)

    return inferred_sequences


def infer_sequences_without_context(bam_file, genome_dict, min_perc_identity, max_internal_softclip_prop):
    bam = pysam.AlignmentFile(bam_file, 'rb')
    keep_reads = prefilter_nocontext_reads(bam, genome_dict, min_perc_identity)
    keep_pairs = get_pairs(keep_reads)
    keep_pairs = filter_nocontext_pairs(keep_pairs, max_internal_softclip_prop)

    for pair_id in keep_pairs:
        keep_pairs[pair_id] = keep_best_alignment_score(keep_pairs[pair_id])

    inferred_sequences = defaultdict(list)
    for pair_id in keep_pairs:
        inferred_sequences[pair_id] = get_inferred_sequences(keep_pairs[pair_id], genome_dict, add_softclipped_bases=True)

    return inferred_sequences


def prefilter_context_reads(bam, genome_dict, min_perc_identity):
    keep_reads = defaultdict(lambda: defaultdict(list))

    for read in bam:
        if read.query_name.count('_') == 1:
            continue

        if pysamtools.get_perc_identity(read) < min_perc_identity:
            continue

        if not read.is_reverse:
            if sctools.is_left_softclipped_strict(read):
                continue

        if read.is_reverse:
            if sctools.is_right_softclipped_strict(read):
                continue

        pair_id, width, flank_id = read.query_name.split('_')

        keep_reads[pair_id+'_'+width][read.reference_name].append(read)

    return keep_reads


def prefilter_nocontext_reads(bam, genome_dict, min_perc_identity):
    keep_reads = defaultdict(lambda: defaultdict(list))

    for read in bam:
        if read.query_name.count('_') != 1:
            continue
        if pysamtools.get_perc_identity(read) < min_perc_identity:
            continue

        if not read.is_reverse:
            if sctools.is_left_softclipped_strict(read) and \
                sctools.get_left_softclip_length(read) > 1 and \
                sctools.left_softclipped_position(read) >= 0:
                continue

        if read.is_reverse:
            if sctools.is_right_softclipped_strict(read) and \
                sctools.get_right_softclip_length(read) > 1 and \
                sctools.right_softclipped_position(read) < len(genome_dict[read.reference_name]):
                continue

        pair_id, flank_id = read.query_name.split('_')

        keep_reads[pair_id][read.reference_name].append(read)

    return keep_reads


def get_pairs(reads):

    keep_pairs = defaultdict(list)
    for pair in reads:
        for ref in reads[pair]:
            sorted_reads = sorted(reads[pair][ref], key=get_start)

            current_forward_read = None
            for read in sorted_reads:
                if not read.is_reverse:
                    current_forward_read = read
                elif current_forward_read is not None:
                        keep_pairs[pair].append((current_forward_read, read))
                        current_forward_read = None
    return keep_pairs


def filter_nocontext_pairs(pairs, max_internal_softclip_prop):

    keep_pairs = defaultdict(list)
    for pair_id in pairs:
        for read1, read2 in pairs[pair_id]:

            if sctools.is_right_softclipped_strict(read1) and \
                read1.reference_end < read2.reference_end and \
                sctools.right_softclip_proportion(read1) > max_internal_softclip_prop:
                continue

            if sctools.is_left_softclipped_strict(read2) and \
                read2.reference_start > read1.reference_start and \
                sctools.left_softclip_proportion(read2) > max_internal_softclip_prop:
                continue

            keep_pairs[pair_id].append((read1, read2))

    return keep_pairs


def filter_context_pairs(pairs, max_internal_softclip_prop):
    keep_pairs = defaultdict(list)
    for pair_id in pairs:
        for read1, read2 in pairs[pair_id]:

            if sctools.is_right_softclipped_strict(read1) and \
                sctools.right_softclip_proportion(read1) > max_internal_softclip_prop:
                continue

            if sctools.is_left_softclipped_strict(read2) and \
                sctools.left_softclip_proportion(read2) > max_internal_softclip_prop:
                continue

            keep_pairs[pair_id].append((read1, read2))

    return keep_pairs


def get_start(read):
    if read.is_reverse:
        return read.reference_end
    else:
        return read.reference_start


def keep_best_alignment_score(reads):
    keep_reads = []

    best_score = 0
    for read1, read2 in reads:
        read_combined_alignment_score = read1.get_tag('AS') + read2.get_tag('AS')
        if read_combined_alignment_score > best_score:
            best_score = read_combined_alignment_score

    for read1, read2 in reads:
        read_combined_alignment_score = read1.get_tag('AS') + read2.get_tag('AS')
        if read_combined_alignment_score == best_score:
            keep_reads.append((read1, read2))

    return keep_reads


def write_flanks_to_align_to_assembly(pairs, bam, genome_dict, tmp_dir):
    fasta_prefix = join(tmp_dir, 'mustache.inferseq_assembly.' + str(randint(0, 1e20)))

    flanks = get_flanks(pairs)
    context_flanks = get_flanks_in_context(pairs, bam, genome_dict)

    writeflanks = flanks+context_flanks
    fastatools.write_flanks_to_unpaired_fasta(writeflanks, fasta_prefix)

    return fasta_prefix


def write_flanks_to_align_to_reference(pairs, tmp_dir):
    fasta_prefix = join(tmp_dir, 'mustache.inferseq_assembly.' + str(randint(0, 1e20)))

    writeflanks = get_flanks(pairs)

    fastatools.write_flanks_to_unpaired_fasta(writeflanks, fasta_prefix)

    return fasta_prefix


def get_flanks(pairs):
    flanks = []
    for index, p in pairs.iterrows():
        pair_id, seq_5p, seq_3p = p['pair_id'], p['seq_5p'], p['seq_3p']
        seq_5p, seq_3p = seq_5p.upper().rstrip('N'), seq_3p.upper().lstrip('N')
        flanks.append({'pair_id':str(pair_id),
                       'seq_5p': seq_5p,
                       'seq_3p': seq_3p})
    return flanks

def get_flanks_in_context(pairs, bam, genome_dict, context_width=25):
    context_flanks = []
    for index, p in pairs.iterrows():
        pair_id, contig_name, pos_5p, pos_3p, seq_5p, seq_3p = p['pair_id'], p['contig'], p['pos_5p'], p['pos_3p'], p['seq_5p'], p['seq_3p']
        seq_5p, seq_3p = seq_5p.rstrip('N'), seq_3p.lstrip('N')
        context_5p = get_sequence_context(bam, genome_dict[contig_name], contig_name, pos_5p-context_width, pos_5p)
        context_3p = get_sequence_context(bam, genome_dict[contig_name], contig_name, pos_3p+1, pos_3p+1+context_width)

        context_flanks.append({'pair_id': str(pair_id)+'_'+str(context_width),
                       'seq_5p': context_5p+seq_5p,
                       'seq_3p': seq_3p+context_3p})

    return context_flanks

def get_inferred_sequences(pairs, genome_dict, add_softclipped_bases=False):

    inferred_sequences = []
    for read1, read2 in pairs:
        if read1.query_name.count('_') == 2:
            context_width = int(read1.query_name.split('_')[-2])
            name = read1.reference_name + ':' + str(read1.reference_start+context_width) + '-' + str(read2.reference_end-context_width)

            inferred_sequence = genome_dict[read1.reference_name][read1.reference_start+context_width:read2.reference_end-context_width]

            if add_softclipped_bases:
                inferred_sequence = sctools.left_softclipped_sequence_strict(read1) + inferred_sequence + sctools.right_softclipped_sequence_strict(read2)

            if read1.query_name.split('_')[-1] == '2':
                inferred_sequence = misc.revcomp(inferred_sequence)

        else:
            name = read1.reference_name + ':' + str(read1.reference_start) + '-' + str(read2.reference_end)
            inferred_sequence = genome_dict[read1.reference_name][read1.reference_start:read2.reference_end]

            if add_softclipped_bases:
                inferred_sequence = sctools.left_softclipped_sequence_strict(read1) + inferred_sequence + sctools.right_softclipped_sequence_strict(read2)

            if read1.query_name.split('_')[-1] == '2':
                inferred_sequence = misc.revcomp(inferred_sequence)

        inferred_sequences.append((name, len(inferred_sequence), inferred_sequence))

    return inferred_sequences

def get_sequence_context(bam, contig, contig_name, start, end):

    add_start_n = 0
    add_end_n = 0
    if start < 0:
        add_start_n = abs(start)
        start = 0
    if end > len(contig):
        add_end_n = end - len(contig)
        end = len(contig)

    reference_sequence = 'N' * add_start_n + contig[start:end] + 'N' * add_end_n
    sequence_context_dict = initialize_sequence_context(reference_sequence, start, end)

    for read in bam.fetch(contig_name, start, end):

        ref_positions = read.get_reference_positions(full_length=True)
        read_query = read.query_sequence
        read_qualities = read.query_qualities

        i = 0
        for pos in ref_positions:
            if pos is not None and pos in sequence_context_dict:
                sequence_context_dict[pos][read_query[i]] += read_qualities[i]
            i += 1

    consensus = get_consensus_context(sequence_context_dict)
    return consensus

def initialize_sequence_context(target_region, expanded_start, expanded_end):
    target_region_reads = defaultdict(lambda: defaultdict(int))
    i = 0
    for pos in range(expanded_start, expanded_end):
        target_region_reads[pos][target_region[i]] += 1
        i += 1
    return target_region_reads

def get_consensus_context(sequence_context_dict):
    start, end = min(sequence_context_dict.keys()), max(sequence_context_dict.keys())+1
    consensus = ''
    for pos in range(start, end):
        best_qual = 0
        best_nuc = ''
        for nuc in sequence_context_dict[pos]:
            qual = sequence_context_dict[pos][nuc]
            if qual > best_qual:
                best_qual = qual
                best_nuc = nuc
        consensus += best_nuc
    return consensus

def index_genome(inferseq_assembly):
    if not bowtie2tools.genome_is_indexed(inferseq_assembly):
        logger.info("Indexing inferseq assembly...")
        bowtie2tools.index_genome(inferseq_assembly)
    logger.info("Genome has been indexed...")


def handle_empty_pairsfile(pairs, output_file):
    if pairs.shape[0] == 0:
        outfile = pd.DataFrame(columns=['pair_id', 'method', 'loc', 'inferred_seq_length', 'inferred_seq'])

        if not output_file:
            output_file = 'mustache.inferseq_assembly.tsv'

        outfile.to_csv(output_file, sep='\t', index=False)
        logger.info("Empty pairs file, exiting...")
        sys.exit()
