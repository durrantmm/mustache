import sys
from multiprocessing import Pool
import misc
import sys
from collections import defaultdict


class KmerNode:
    next_nodes = None
    count = None
    read_names = None

    def __init__(self):
        self.next_nodes = set()
        self.count = 0
        self.read_names = set()


class KmerGraph:
    graph = defaultdict(lambda: KmerNode())

    def __init__(self):
        pass

    def add_kmer(self, kmer, read_name, next_kmer=None):
        self.graph[kmer].count += 1
        self.graph[kmer].read_names.add(read_name)

        if next_kmer:
            self.graph[kmer].next_nodes.add(next_kmer)

    def __getitem__(self, x):
        return self.graph[x]

    def __str__(self):

        out = ''
        for kmer in self.graph:
            count = self.graph[kmer].count
            next_nodes = ','.join(list(self.graph[kmer].next_nodes))[0:45] + '...'
            read_names = ','.join(list(self.graph[kmer].read_names))[0:45] + '...'

            next_nodes += ' ' * (48 - len(next_nodes))
            read_names += ' ' * (48 - len(read_names))
            out += '\t'.join(map(str, [kmer, count, next_nodes, read_names])) + '\n'

        return out


class InsertionSequenceAssembler:

    starting_node = None
    kmer_graph = KmerGraph()
    softclipped_reads = None
    unmapped_reads = None
    total_kmer_count = None

    def __init__(self, softclipped_reads, unmapped_reads, kmer_size=31):

        self.softclipped_reads = softclipped_reads
        self.unmapped_reads = unmapped_reads
        self.total_kmer_count = defaultdict(int)

        first_kmer = defaultdict(int)

        for r in softclipped_reads:
            read_name = r[0]
            read_seq = r[1]

            if len(read_seq) <= kmer_size:
                continue

            kmers = list(self.kmer_composition(read_seq, kmer_size))
            first_kmer[kmers[0]] += 1

            for i in range(len(kmers)):
                self.total_kmer_count[read_name] += 1
                try:
                    self.kmer_graph.add_kmer(kmers[i], read_name, kmers[i + 1])

                except IndexError:
                    self.kmer_graph.add_kmer(kmers[i], read_name)

        # Get starting kmer
        self.starting_node = max(list(first_kmer.items()), key=lambda x: x[1])[0]

        # Now add all of the other kmers.
        for r in self.unmapped_reads:
            read_name = r[0]
            forward_read_seq = r[1]
            reverse_read_seq = misc.revcomp(r[1])

            if len(forward_read_seq) <= kmer_size:
                continue

            forward_kmers = list(self.kmer_composition(forward_read_seq, kmer_size))
            reverse_kmers = list(self.kmer_composition(reverse_read_seq, kmer_size))

            for i in range(len(forward_kmers)):
                self.total_kmer_count[read_name] += 1
                try:
                    self.kmer_graph.add_kmer(forward_kmers[i], read_name, forward_kmers[i + 1])
                    self.kmer_graph.add_kmer(reverse_kmers[i], read_name, reverse_kmers[i + 1])

                except IndexError:
                    self.kmer_graph.add_kmer(forward_kmers[i], read_name)
                    self.kmer_graph.add_kmer(reverse_kmers[i], read_name)


    def assemble(self, trim_by_softclips=True, max_len=3000):
        read_name_count = defaultdict(int)

        current_node = self.starting_node

        assembly = current_node
        while len(self.kmer_graph[current_node].next_nodes) > 0:
            for name in self.kmer_graph[current_node].read_names:
                read_name_count[name] += 1
            next_node_count = 0
            next_node = None

            for node in self.kmer_graph[current_node].next_nodes:
                if self.kmer_graph[node].count > next_node_count:
                    next_node = node
                    next_node_count = self.kmer_graph[node].count

            assembly += next_node[-1]
            current_node = next_node

            if len(assembly) >= max_len:
                print("REACHED MAXIMUM LENGTH")
                break

        if trim_by_softclips:
            start_sites = defaultdict(int)
            for r in self.softclipped_reads:
                sc_seq = r[2]
                best_score, best_score_mismatches, best_r_start, best_r_end = self.get_best_sliding_alignment(sc_seq,
                                                                                                              assembly)
                start_sites[best_r_start] += 1
            start_site = max(list(start_sites.items()), key=lambda x: x[1])[0]
            assembly = assembly[start_site:]

        return assembly


    def kmer_composition(self, read, kmer_size):
        for i in range(len(read) - kmer_size):
            yield read[i:(i + kmer_size)]


    def get_best_sliding_alignment(self, query_read, ref_read):

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

        return best_score, best_score_mismatches, best_r_start, best_r_end


def get_simple_alignment_score(read1, read2):
    score = 0
    for i in range(min([len(read1), len(read2)])):
        if read1[i] == read2[i]:
            score += 1
        else:
            score -= 1

    return (score)


def get_max_sliding_alignment_score(query_read, ref_read):
    max_len = min(len(query_read), len(ref_read))
    max_len_query, max_len_ref = False, False

    if max_len == len(query_read):
        max_len_query = True
    else:
        max_len_ref = True

    best_score = -sys.maxsize
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

    return best_score, best_q_end - best_q_start


def assemble_flank(softclipped, unmapped, kmer_sizes=[21,31,41,51,61]):

    longest_assembly_length = 0
    longest_assembly = ''
    for kmer_size in kmer_sizes:
        assembler = InsertionSequenceAssembler(softclipped, unmapped, kmer_size)
        assembly = assembler.assemble()
        if len(assembly) > longest_assembly_length:
            longest_assembly_length = len(assembly)
            longest_assembly = assembly
    return longest_assembly




