import sys
from mustache import misc
import click
from Bio.pairwise2 import format_alignment
from Bio import pairwise2

class QuickMapper:

    reference_sequence = None
    kmer_size = None
    kmer_dict = None

    def __init__(self, reference_sequence, kmer_size=10, min_alignment_score=11):
        self.reference_sequence = reference_sequence
        self.kmer_size = kmer_size
        self.min_alignment_score = min_alignment_score
        self.kmer_dict = dict()

        kmers = list(self.kmer_composition(reference_sequence))

        for pos in range(len(kmers)):

            kmer = kmers[pos]

            if kmer not in self.kmer_dict.keys():
                self.kmer_dict[kmer] = dict()
                self.kmer_dict[kmer]['count'] = 1
                self.kmer_dict[kmer]['positions'] = [pos]
            else:
                self.kmer_dict[kmer]['count'] += 1
                if pos not in self.kmer_dict[kmer]['positions']:
                    self.kmer_dict[kmer]['positions'].append(pos)

    def print_kmer_dict(self):
        for key in self.kmer_dict:
            print(key, self.kmer_dict[key])

    def kmer_composition(self, read):
        for i in range(len(read) - self.kmer_size):
            yield read[i:(i + self.kmer_size)]

    def maps_to_reference(self, seq):
        forward_seq = seq
        reverse_seq = misc.revcomp(seq)

        kmers_match = False
        best_mapping_sequence = None
        best_mapping_positions = None
        total_matching_kmers = 0

        for sequence in [forward_seq, reverse_seq]:
            num_matching_kmers = 0
            mapping_positions = []

            kmers = list(self.kmer_composition(sequence))

            for kmer in kmers:
                if kmer in self.kmer_dict:
                    mapping_positions.append(self.kmer_dict[kmer]['positions'])
                    num_matching_kmers += 1
                else:
                    mapping_positions.append([])

            if kmers_match:
                if num_matching_kmers > total_matching_kmers:
                    best_mapping_sequence = sequence
                    best_mapping_positions = mapping_positions
                    total_matching_kmers = num_matching_kmers
            elif num_matching_kmers > 0:
                best_mapping_sequence = sequence
                best_mapping_positions = mapping_positions
                total_matching_kmers = num_matching_kmers
                kmers_match = True

        if kmers_match:

            sequence_start, sequence_end = None, None
            reference_start, reference_end = None, None

            for i in range(len(best_mapping_positions)):
                if len(best_mapping_positions[i]) == 1:
                    sequence_start = i
                    reference_start = best_mapping_positions[i][0]
                    break
            for i in reversed(range(len(best_mapping_positions))):
                if len(best_mapping_positions[i]) == 1:
                    sequence_end = i+1
                    reference_end = best_mapping_positions[i][0] + 1
                    break

            clipped_sequence = best_mapping_sequence[sequence_start:sequence_end]
            clipped_reference = self.reference_sequence[reference_start:reference_end]

            left_score = left_alignment_score(clipped_sequence, clipped_reference)
            right_score = right_alignment_score(clipped_sequence, clipped_reference)

            if max([left_score, right_score]) >= self.min_alignment_score:
                return True

        return False

    def maps_within_reference(self, seq):
        forward_seq = seq
        reverse_seq = misc.revcomp(seq)

        kmers_match = False
        best_mapping_sequence = None
        best_mapping_positions = None
        total_matching_kmers = -1

        for sequence in [forward_seq, reverse_seq]:
            num_matching_kmers = 0
            mapping_positions = []

            kmers = list(self.kmer_composition(sequence))

            for kmer in kmers:
                if kmer in self.kmer_dict:
                    mapping_positions.append(self.kmer_dict[kmer]['positions'])
                    num_matching_kmers += 1
                else:
                    mapping_positions.append([])

            if kmers_match:
                if num_matching_kmers > total_matching_kmers:
                    best_mapping_sequence = sequence
                    best_mapping_positions = mapping_positions
                    total_matching_kmers = num_matching_kmers

            elif num_matching_kmers > 0:
                best_mapping_sequence = sequence
                best_mapping_positions = mapping_positions
                total_matching_kmers = num_matching_kmers
                kmers_match = True

        if kmers_match:

            sequence_start, sequence_end = None, None
            reference_start, reference_end = None, None

            for i in range(len(best_mapping_positions)):
                if len(best_mapping_positions[i]) == 1:
                    sequence_start = i
                    reference_start = best_mapping_positions[i][0]
                    break

            for i in reversed(range(len(best_mapping_positions))):
                if len(best_mapping_positions[i]) == 1:
                    sequence_end = i+1+self.kmer_size
                    reference_end = best_mapping_positions[i][0] + 1 + self.kmer_size
                    break

            if sequence_start is None or sequence_end is None or reference_start is None or reference_end is None:
                return False

            # Now check that the sequence is well contained within the sequencing reference inseq.
            if sequence_start > 0 and reference_start == 0:
                return False

            if len(best_mapping_sequence) - sequence_end > 0 and len(self.reference_sequence) - reference_end == 0:
                return False

            # Now continue on to check the score of the alignment
            clipped_sequence = best_mapping_sequence[sequence_start:sequence_end]
            clipped_reference = self.reference_sequence[reference_start:reference_end]

            left_score = left_alignment_score(clipped_sequence, clipped_reference)
            right_score = right_alignment_score(clipped_sequence, clipped_reference)

            if max([left_score, right_score]) >= self.min_alignment_score:
                return True

        return False

def left_alignment_score(read1, read2):
    score = 0
    for i in range(min([len(read1), len(read2)])):
        if read1[i] == read2[i]:
            score += 1
        else:
            score -= 1

    return (score)


def right_alignment_score(read1, read2):
    score = 0

    read1 = read1[::-1]
    read2 = read2[::-1]

    for i in range(min([len(read1), len(read2)])):
        if read1[i] == read2[i]:
            score += 1
        else:
            score -= 1

    return (score)


def get_best_sliding_alignment(query_read, ref_read):
    max_len = min(len(query_read), len(ref_read))
    max_len_query, max_len_ref = False, False

    if max_len == len(query_read):
        max_len_query = True
    else:
        max_len_ref = True

    best_score = -sys.maxsize
    best_score_mismatches = 0
    best_r_start = 0
    best_r_end = 0
    best_q_start = 0
    best_q_end = 0

    q_start = len(query_read) - 1
    q_end = len(query_read)
    r_start = 0
    r_end = 1

    reached_max = False

    while q_start != q_end:
        score = 0
        mismatches = 0
        for i in range(q_end - q_start):
            if query_read[q_start:q_end][i] == ref_read[r_start:r_end][i]:
                score += 1
            else:
                score -= 1
                mismatches += 1

        if score > best_score:
            best_score = score
            best_score_mismatches = mismatches
            best_r_start = r_start
            best_r_end = r_end
            best_q_start = q_start
            best_q_end = q_end

        if q_end - q_start == max_len:
            reached_max = True

        if not reached_max:
            q_start -= 1
            r_end += 1

        elif reached_max and max_len_query:
            if r_end == len(ref_read):
                r_start += 1
                q_end -= 1
            else:
                r_start += 1
                r_end += 1

        elif reached_max and max_len_ref:
            if q_start == 0:
                r_start += 1
                q_end -= 1
            else:
                q_start -= 1
                q_end -= 1

    return best_score, best_score_mismatches, best_r_start, best_r_end, best_q_start, best_q_end


def merge_flank_assemblies(right_inseq_assembly, left_inseq_assembly, min_overlap_score = 10, mismatch_prop= 1.0 / 10.0, verbose=True):
    align_info = get_best_sliding_alignment(right_inseq_assembly, left_inseq_assembly)
    best_score, best_score_mismatches, best_r_start, best_r_end, best_q_start, best_q_end = align_info
    align_length = best_r_end - best_r_start

    merged_assembly = None

    if best_score >= min_overlap_score and best_score_mismatches / float(align_length) < mismatch_prop:
        if verbose:
            click.echo("They merged!")
            alignments = pairwise2.align.localms(right_inseq_assembly, left_inseq_assembly, 1, -1, -5, -5)
            print(format_alignment(*alignments[0]))


        if best_q_start > 0 and best_q_end == len(right_inseq_assembly) and \
            best_r_start == 0 and best_r_end != len(left_inseq_assembly):
            # XXXXXXXXXXXX
            #        XXXXXXXXXXXX
            # XXXXXXXXXXXXXXXXXXX
            merged_assembly = right_inseq_assembly + left_inseq_assembly[best_r_end:]

        elif best_q_start == 0 and best_q_end == len(right_inseq_assembly) and \
            best_r_start == 0 and best_r_end == len(left_inseq_assembly):
            # XXXXXXXXXXXXXXX
            # XXXXXXXXXXXXXXX
            # XXXXXXXXXXXXXXX
            merged_assembly = right_inseq_assembly

        elif best_q_start == 0 and best_q_end == len(right_inseq_assembly) and \
            best_r_start > 0 and best_r_end != len(left_inseq_assembly):
            #       XXXXXXXXX
            # XXXXXXXXXXXXXXXXXXX
            #       XXXXXXXXXXXXX
            merged_assembly = left_inseq_assembly[best_r_start:]

        elif best_q_start > 0 and best_q_end != len(right_inseq_assembly) and \
            best_r_start == 0 and best_r_end == len(left_inseq_assembly):
            # XXXXXXXXXXXXXXXXXXX
            #      XXXXXXXXX
            # XXXXXXXXXXXXXX
            merged_assembly = right_inseq_assembly[:best_q_end]

        elif best_q_start == 0 and best_q_end != len(right_inseq_assembly) and \
            best_r_start > 0 and best_r_end == len(left_inseq_assembly):
            #        XXXXXXXXXXXX
            # XXXXXXXXXXXX
            #        XXXXX
            merged_assembly = right_inseq_assembly[best_q_start:best_q_end]
        elif best_q_start == 0 and best_q_end == len(right_inseq_assembly) and \
                best_r_start > 0 and best_r_end == len(left_inseq_assembly):
            #     XXXXXXXXXXXXXXX
            # XXXXXXXXXXXXXXXXXXX
            #     XXXXXXXXXXXXXXX
            merged_assembly = right_inseq_assembly
        elif best_q_start == 0 and best_q_end < len(right_inseq_assembly) and \
            best_r_start == 0 and best_r_end == len(left_inseq_assembly):
            # XXXXXXXXXXXXXXXXXXX
            # XXXXXXXXXXXXXXX
            # XXXXXXXXXXXXXXX
            merged_assembly = left_inseq_assembly
        elif best_q_start == 0 and best_q_end == len(right_inseq_assembly) and \
            best_r_start == 0 and best_r_end < len(left_inseq_assembly):
            # XXXXXXXXXXXXXXX
            # XXXXXXXXXXXXXXXXXXX
            # XXXXXXXXXXXXXXXXXXX
            merged_assembly = left_inseq_assembly
        else:
            print(right_inseq_assembly)
            print(left_inseq_assembly)

            print('Length Right:', len(right_inseq_assembly))
            print('Length Left:', len(left_inseq_assembly))
            print(best_q_start, best_q_end, best_r_start, best_r_end)
            print("I MISSED SOMETHING")
            sys.exit()

    return right_inseq_assembly, left_inseq_assembly, merged_assembly