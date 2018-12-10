import sys
from mustache import pysamtools, sctools, misc, flanktrie
import pysam
from collections import defaultdict
import pandas as pd
from jellyfish import levenshtein_distance
import itertools
from scipy.sparse.csgraph import connected_components
import numpy as np

import pygogo as gogo
verbose = True
logger = gogo.Gogo(__name__, verbose=False).logger



def _findflanks(bamfile, min_softclip_length, min_softclip_count, min_alignment_quality, min_alignment_inner_length,
                min_distance_to_mate, min_softclip_ratio, max_indel_ratio, min_count_consensus, output_file=None):

    bam = pysam.AlignmentFile(bamfile, 'rb')

    softclip_parser = SoftclipParser(
        bam, verbose=True,
        min_softclip_length=min_softclip_length,
        min_softclip_count=min_softclip_count,
        min_alignment_quality=min_alignment_quality,
        min_alignment_inner_length=min_alignment_inner_length,
        min_distance_to_mate=min_distance_to_mate,
        min_softclip_ratio=min_softclip_ratio,
        max_indel_ratio=max_indel_ratio,
        min_count_consensus=min_count_consensus
    )

    softclip_parser.parse_softclips()
    softclip_parser.filter_softclips_minlength()
    softclip_parser.filter_softclips_mincount()
    softclip_parser.filter_softclips_mindistance()
    softclip_parser.parse_unclipped_read_info()
    softclip_parser.filter_softclips_count_ratios()
    softclip_parser.filter_softclips_mindistance()
    softclip_parser.make_consensus_sequences()
    softclip_parser.filter_consensus_sequences_minlength()
    softclip_parser.filter_consensus_sequences_mincount()
    softclip_parser.filter_multiple_consensus_sequences()
    softclip_parser.filter_softclips_mindistance()

    final_df = softclip_parser.make_dataframe()
    final_df.index.names = ['flank_id']
    final_df.index = final_df.index + 1

    if output_file:
        logger.info("Saving results to file %s" % output_file)
        final_df.to_csv(output_file, sep='\t')

    return final_df


class SoftclipParser:

    softclipped_sites = None
    bam = None
    contig_lengths = None
    verbose = None

    min_alignment_quality = None
    min_alignment_inner_length = None
    min_softclip_length = None
    min_softclip_count = None
    min_distance_to_mate = None
    min_softclip_ratio = None
    max_indel_ratio = None
    min_count_consensus = None

    def __init__(self, bam, verbose=True, min_alignment_quality=20, min_alignment_inner_length=21,
                 min_softclip_length=4, min_softclip_count=4, min_distance_to_mate=22,
                 min_softclip_ratio=0.15, max_indel_ratio=0.03, min_count_consensus=2):
        self.verbose = verbose
        self.bam = bam
        self.contig_lengths = pysamtools.get_bam_contig_dict(bam)

        self.softclipped_sites = defaultdict(lambda: defaultdict(SoftclipSite))
        self.min_alignment_quality = min_alignment_quality
        self.min_alignment_inner_length = min_alignment_inner_length
        self.min_softclip_length = min_softclip_length
        self.min_softclip_count = min_softclip_count
        self.min_distance_to_mate = min_distance_to_mate
        self.min_softclip_ratio = min_softclip_ratio
        self.max_indel_ratio = max_indel_ratio
        self.min_count_consensus = min_count_consensus

    def parse_softclips(self):

        if self.verbose:
            logger.info("Parsing softclipped sites from provided BAM file...")

        read_count = 0
        for read in self.bam:

            read_count += 1

            if self.verbose and read_count % 100000 == 0:
                logger.info("\tAfter checking %d reads, %d softclipped sites found..." % (read_count, self.count_softclips()))
                pass

            if not self.passes_read_filters(read):
                continue

            if sctools.is_right_softclipped_lenient(read):
                self.parse_right_softclipped_read(read)

            if sctools.is_left_softclipped_lenient(read):
                self.parse_left_softclipped_read(read)


    def filter_softclips_minlength(self):
        filtered_softclipped_sites = defaultdict(lambda: defaultdict(SoftclipSite))

        for contig in self.softclipped_sites:
            for pos in self.softclipped_sites[contig]:

                if not self.softclipped_sites[contig][pos].meets_minlength_5p:
                    self.softclipped_sites[contig][pos].keep_softclips_5p = False

                if not self.softclipped_sites[contig][pos].meets_minlength_3p:
                    self.softclipped_sites[contig][pos].keep_softclips_3p = False

                if self.softclipped_sites[contig][pos].keep_softclips_5p or \
                    self.softclipped_sites[contig][pos].keep_softclips_3p:
                    filtered_softclipped_sites[contig][pos] = self.softclipped_sites[contig][pos]

        self.softclipped_sites = filtered_softclipped_sites

        if verbose:
            logger.info("After filtering by minimum softclip length of %d, %d sites remain" % (self.min_softclip_length, self.count_softclips()))
            pass


    def filter_softclips_mincount(self):
        filtered_softclipped_sites = defaultdict(lambda: defaultdict(SoftclipSite))

        for contig in self.softclipped_sites:
            for pos in self.softclipped_sites[contig]:
                if self.softclipped_sites[contig][pos].get_softclip_5p_count() < self.min_softclip_count:
                    self.softclipped_sites[contig][pos].keep_softclips_5p = False

                if self.softclipped_sites[contig][pos].get_softclip_3p_count() < self.min_softclip_count:
                    self.softclipped_sites[contig][pos].keep_softclips_3p = False

                if self.softclipped_sites[contig][pos].keep_softclips_5p or \
                    self.softclipped_sites[contig][pos].keep_softclips_3p:
                    filtered_softclipped_sites[contig][pos] = self.softclipped_sites[contig][pos]

        self.softclipped_sites = filtered_softclipped_sites

        if verbose:
            logger.info("After filtering by minimum softclipped read count of %d, %d sites remain" % (
                self.min_softclip_count, self.count_softclips()))
            pass

    def filter_softclips_mindistance(self):

        softclip_5p_positions = self.get_softclip_5p_positions()
        softclip_3p_positions = self.get_softclip_3p_positions()

        filtered_softclipped_sites = defaultdict(lambda: defaultdict(SoftclipSite))

        for contig in self.softclipped_sites:
            for pos in self.softclipped_sites[contig]:
                softclip_site = self.softclipped_sites[contig][pos]

                if softclip_site.get_softclip_5p_count() > 0:

                    if not self.has_nearby_3p_mate(pos, softclip_3p_positions[contig]):
                        softclip_site.keep_softclips_5p = False

                if softclip_site.get_softclip_3p_count() > 0:

                    if not self.has_nearby_5p_mate(pos, softclip_5p_positions[contig]):
                        softclip_site.keep_softclips_3p = False

                if self.softclipped_sites[contig][pos].keep_softclips_5p or \
                    self.softclipped_sites[contig][pos].keep_softclips_3p:
                    filtered_softclipped_sites[contig][pos] = softclip_site

        self.softclipped_sites = filtered_softclipped_sites

        if verbose:
            logger.info("After filtering by minimum nearest mate distance %d, %d sites remain" % (
                self.min_distance_to_mate, self.count_softclips()))
            pass


    def filter_softclips_count_ratios(self):

        filtered_softclipped_sites = defaultdict(lambda: defaultdict(SoftclipSite))

        for contig in self.softclipped_sites:
            for pos in self.softclipped_sites[contig]:
                softclip_site = self.softclipped_sites[contig][pos]

                if softclip_site.keep_softclips_5p:

                    if softclip_site.get_softclip_ratio_5p() < self.min_softclip_ratio or \
                        softclip_site.get_indel_ratio_5p() > self.max_indel_ratio:
                            softclip_site.keep_softclips_5p = False

                if softclip_site.keep_softclips_3p:

                    if softclip_site.get_softclip_ratio_3p() < self.min_softclip_ratio or \
                        softclip_site.get_indel_ratio_3p() > self.max_indel_ratio:
                        softclip_site.keep_softclips_3p = False

                if softclip_site.keep_softclips_5p or softclip_site.keep_softclips_3p:

                    if softclip_site.get_downstream_deletion_ratio() > self.max_indel_ratio or \
                       softclip_site.get_upstream_deletion_ratio() > self.max_indel_ratio:
                        softclip_site.keep_softclips_5p = False
                        softclip_site.keep_softclips_3p = False

                if self.softclipped_sites[contig][pos].keep_softclips_5p or \
                    self.softclipped_sites[contig][pos].keep_softclips_3p:
                    filtered_softclipped_sites[contig][pos] = softclip_site

        self.softclipped_sites = filtered_softclipped_sites

        if verbose:
            logger.info("After filtering by minimum softclip ratio of %f and a "
                        "maximum indel ratio of %f, %d sites remain" % (
                self.min_softclip_ratio, self.max_indel_ratio, self.count_softclips()))
            pass

    def filter_consensus_sequences_minlength(self):

        filtered_softclipped_sites = defaultdict(lambda: defaultdict(SoftclipSite))

        for contig in self.softclipped_sites:
            for pos in self.softclipped_sites[contig]:
                softclip_site = self.softclipped_sites[contig][pos]

                if softclip_site.keep_softclips_5p:

                    consensus_seqs = []
                    for consensus in softclip_site.consensus_sequences_5p:

                        if consensus[0] >= self.min_softclip_length:
                            consensus_seqs.append(consensus)

                    softclip_site.consensus_sequences_5p = consensus_seqs
                    if len(consensus_seqs) == 0:
                        softclip_site.keep_softclips_5p = False

                if softclip_site.keep_softclips_3p:

                    consensus_seqs = []
                    for consensus in softclip_site.consensus_sequences_3p:

                        if consensus[0] >= self.min_softclip_length:
                            consensus_seqs.append(consensus)

                    softclip_site.consensus_sequences_3p = consensus_seqs
                    if len(consensus_seqs) == 0:
                        softclip_site.keep_softclips_3p = False

                if self.softclipped_sites[contig][pos].keep_softclips_5p or \
                        self.softclipped_sites[contig][pos].keep_softclips_3p:
                    filtered_softclipped_sites[contig][pos] = softclip_site

        self.softclipped_sites = filtered_softclipped_sites

        if verbose:
            logger.info("After filtering consensus sequences by a minimum length of %d, %d flank sequences remain" % (
                            self.min_softclip_count, self.count_consensus_seqs())
                        )

    def filter_consensus_sequences_mincount(self):

        filtered_softclipped_sites = defaultdict(lambda: defaultdict(SoftclipSite))

        for contig in self.softclipped_sites:
            for pos in self.softclipped_sites[contig]:
                softclip_site = self.softclipped_sites[contig][pos]

                if softclip_site.keep_softclips_5p:

                    consensus_seqs = []
                    for consensus in softclip_site.consensus_sequences_5p:

                        if consensus[0] >= self.min_softclip_count:
                            consensus_seqs.append(consensus)
                    softclip_site.consensus_sequences_5p = consensus_seqs
                    if len(consensus_seqs) == 0:
                        softclip_site.keep_softclips_5p = False

                if softclip_site.keep_softclips_3p:

                    consensus_seqs = []
                    for consensus in softclip_site.consensus_sequences_3p:

                        if consensus[0] >= self.min_softclip_count:
                            consensus_seqs.append(consensus)
                    softclip_site.consensus_sequences_3p = consensus_seqs
                    if len(consensus_seqs) == 0:
                        softclip_site.keep_softclips_3p = False

                if self.softclipped_sites[contig][pos].keep_softclips_5p or \
                        self.softclipped_sites[contig][pos].keep_softclips_3p:
                    filtered_softclipped_sites[contig][pos] = softclip_site

        self.softclipped_sites = filtered_softclipped_sites

        if verbose:
            logger.info("After filtering consensus sequences by a minimum softclip count of %d, %d flank sequences remain" % (
                self.min_softclip_count, self.count_consensus_seqs())
                        )

    def filter_multiple_consensus_sequences(self):

        filtered_softclipped_sites = defaultdict(lambda: defaultdict(SoftclipSite))

        for contig in self.softclipped_sites:
            for pos in self.softclipped_sites[contig]:
                softclip_site = self.softclipped_sites[contig][pos]


                if len(softclip_site.consensus_sequences_5p) > 1:
                    softclip_site.keep_softclips_5p = False

                if len(softclip_site.consensus_sequences_3p) > 1:
                    softclip_site.keep_softclips_3p = False

                if self.softclipped_sites[contig][pos].keep_softclips_5p or \
                        self.softclipped_sites[contig][pos].keep_softclips_3p:
                    filtered_softclipped_sites[contig][pos] = softclip_site

        self.softclipped_sites = filtered_softclipped_sites

        if verbose:
            logger.info("After filtering out sites with multiple consensus sequences, %d sites remain." % (
                self.count_softclips())
                        )


    def parse_unclipped_read_info(self):
        if verbose:
            logger.info("Getting unclipped read information near softclipped sites...")
            pass

        for contig in self.softclipped_sites:
            for pos in self.softclipped_sites[contig]:

                runthrough_count, insertion_5p_count, insertion_3p_count, deletion_count = \
                    self.get_unclipped_read_info_at_site(contig, pos)

                self.softclipped_sites[contig][pos].add_runthrough_count(runthrough_count)
                self.softclipped_sites[contig][pos].add_insertion_5p_count(insertion_5p_count)
                self.softclipped_sites[contig][pos].add_insertion_3p_count(insertion_3p_count)
                self.softclipped_sites[contig][pos].add_deletion_count(deletion_count)

                upstream_runthrough_count, upstream_insertion_5p_count, upstream_insertion_3p_count, upstream_deletion_count= None, None, None, None
                downstream_runthrough_count, downstream_insertion_5p_count, downstream_insertion_3p_count, downstream_deletion_count = None, None, None, None

                if pos - 1 >= 0:
                    upstream_runthrough_count, upstream_insertion_5p_count, upstream_insertion_3p_count, upstream_deletion_count = \
                        self.get_unclipped_read_info_at_site(contig, pos + 1)

                if pos + 1 < self.contig_lengths[contig]:
                    downstream_runthrough_count, downstream_insertion_5p_count, downstream_insertion_3p_count, downstream_deletion_count = \
                        self.get_unclipped_read_info_at_site(contig, pos + 1)

                if upstream_deletion_count:
                    self.softclipped_sites[contig][pos].add_upstream_deletion_count(upstream_deletion_count)

                if downstream_deletion_count:
                    self.softclipped_sites[contig][pos].add_downstream_deletion_count(downstream_deletion_count)



    def get_unclipped_read_info_at_site(self, contig, pos):

        runthrough_count = 0
        insertion_5p_count = 0
        insertion_3p_count = 0
        deletion_count = 0

        start = pos-1
        end = pos+2

        if start < 0:
            start = 0
        if end > self.contig_lengths[contig]:
            end = self.contig_lengths[contig]

        for read in self.bam.fetch(contig, start, end):

            if not self.passes_read_filters(read):
                continue

            if sctools.is_softclipped_lenient_at_site(read, contig, pos):
                continue

            processed_read = self.process_aligned_blocks_at_site(pos, read.get_blocks())
            if processed_read is None:
                continue
            elif processed_read == 'runthrough':
                runthrough_count += 1
            elif processed_read == 'insertion_5p':
                insertion_5p_count += 1
            elif processed_read == 'insertion_3p':
                insertion_3p_count += 1
            else:
                deletion_count += 1

        return runthrough_count, insertion_5p_count, insertion_3p_count, deletion_count


    def process_aligned_blocks_at_site(self, position, blocks):

        runthrough = False
        insertion_5p = False
        insertion_3p = False
        deletion = False

        for b in blocks:

            if self.block_overlaps_site(b, position):
                runthrough = True

        if len(blocks) > 1:
            for i in range(len(blocks)-1):
                block1, block2 = blocks[i], blocks[i+1]
                if block1[1] == block2[0]:

                    if position == block1[1]-1:
                        insertion_3p = True
                    elif position == block2[0]:
                        insertion_5p = True
                else:
                    deletion_range = (block1[1], block2[0])
                    if self.block_overlaps_site(deletion_range, position):
                        deletion = True
        if insertion_5p:
            return 'insertion_5p'
        elif insertion_3p:
            return 'insertion_3p'
        elif deletion:
            return 'deletion'
        elif runthrough:
            return 'runthrough'
        else:
            None


    def make_consensus_sequences(self):

        if self.verbose:
            logger.info('Generating consensus sequences from softclipped flanks...')
        for contig in self.softclipped_sites:

            for pos in self.softclipped_sites[contig]:

                site = self.softclipped_sites[contig][pos]

                if site.keep_softclips_5p:

                    reads = site.softclip_5p_reads
                    softclip_consensus = SoftclipConsensus(reads, '5p', pos, self.min_count_consensus)
                    site.consensus_sequences_5p = softclip_consensus.consensus_seqs

                if site.keep_softclips_3p:

                    reads = site.softclip_3p_reads
                    softclip_consensus = SoftclipConsensus(reads, '3p', pos, self.min_count_consensus)
                    site.consensus_sequences_3p = softclip_consensus.consensus_seqs


    def block_overlaps_site(self, block, position):

        if block[0] <= position and block[1] > position:
            return True
        else:
            return False

    def has_nearby_5p_mate(self, pos, positions):
        closest_5p = misc.takeClosestLarger(positions, pos)
        if closest_5p is None:
            return False
        if closest_5p - pos <= self.min_distance_to_mate:
            return True
        else:
            return False

    def has_nearby_3p_mate(self, pos, positions):
        closest_3p = misc.takeClosestSmaller(positions, pos)
        if closest_3p is None:
            return False
        if pos - closest_3p <= self.min_distance_to_mate:
            return True
        else:
            return False

    def get_softclip_5p_positions(self):
        positions = defaultdict(list)

        for contig in self.softclipped_sites:
            for pos in self.softclipped_sites[contig]:
                if self.softclipped_sites[contig][pos].get_softclip_5p_count() > 0:
                    positions[contig].append(pos)

        for contig in positions:
            positions[contig] = sorted(positions[contig])

        return positions

    def get_softclip_3p_positions(self):
        positions = defaultdict(list)

        for contig in self.softclipped_sites:
            for pos in self.softclipped_sites[contig]:
                if self.softclipped_sites[contig][pos].get_softclip_3p_count() > 0:
                    positions[contig].append(pos)

        for contig in positions:
            positions[contig] = sorted(positions[contig])

        return positions

    def passes_read_filters(self, read):
        if read.mapping_quality < self.min_alignment_quality:
            return False
        elif not sctools.read_meets_min_alignment_inner_length(read, self.min_alignment_inner_length):
            return False
        else:
            return True

    def parse_right_softclipped_read(self, read):
        meets_minlength = self.meets_minlength_right(read)
        contig, pos = sctools.right_softclipped_site_lenient(read)
        self.softclipped_sites[contig][pos].add_softclip_5p(read, meets_minlength)


    def parse_left_softclipped_read(self, read):
        meets_minlength = self.meets_minlength_left(read)
        contig, pos = sctools.left_softclipped_site_lenient(read)
        self.softclipped_sites[contig][pos].add_softclip_3p(read, meets_minlength)

    def meets_minlength_right(self, read):
        meets_minlength = False
        if sctools.right_softclip_length(read) >= self.min_softclip_length:
            meets_minlength = True
        return meets_minlength

    def meets_minlength_left(self, read):
        meets_minlength = False
        if sctools.left_softclip_length(read) >= self.min_softclip_length:
            meets_minlength = True
        return meets_minlength


    def count_softclips(self):
        count = 0
        for contig in self.softclipped_sites:
            count += len(self.softclipped_sites[contig])
        return count

    def count_consensus_seqs(self):
        count = 0
        for contig in self.softclipped_sites:
            for pos in self.softclipped_sites[contig]:
                count += len(self.softclipped_sites[contig][pos].consensus_sequences_5p)
                count += len(self.softclipped_sites[contig][pos].consensus_sequences_3p)
        return count

    def make_dataframe(self):

        column_names = ['contig', 'pos', 'orient', 'softclip_count_5p', 'softclip_count_3p', 'runthrough_count',
                        'small_insertion_count_5p', 'small_insertion_count_3p', 'deletion_count',
                        'upstream_deletion_count', 'downstream_deletion_count', 'total_count',
                        'consensus_softclip_count', 'consensus_seq']

        outdata= dict()
        for contig in self.softclipped_sites:
            sorted_positions = sorted(list(self.softclipped_sites[contig].keys()))

            for pos in sorted_positions:
                site = self.softclipped_sites[contig][pos]

                if site.keep_softclips_5p:

                    for consensus_sequence in site.consensus_sequences_5p:

                        outdata[len(outdata)] = [
                            contig, pos, '5p', site.get_softclip_5p_count(), site.get_softclip_3p_count(),
                            site.runthrough_count, site.insertion_5p_count, site.insertion_3p_count,
                            site.deletion_count, site.upstream_deletion_count, site.downstream_deletion_count,
                            site.total_count, consensus_sequence[0], consensus_sequence[1]
                        ]

                if site.keep_softclips_3p:
                    for consensus_sequence in site.consensus_sequences_3p:
                        outdata[len(outdata)] = [
                            contig, pos, '3p', site.get_softclip_5p_count(), site.get_softclip_3p_count(),
                            site.runthrough_count, site.insertion_5p_count, site.insertion_3p_count,
                            site.deletion_count, site.upstream_deletion_count, site.downstream_deletion_count,
                            site.total_count, consensus_sequence[0], consensus_sequence[1]
                        ]

        out_df = pd.DataFrame.from_dict(outdata, orient='index', columns=column_names)

        return out_df


    def print_sites(self):

        print('contig\tpos\ttotal\trunthrough_count\tsoftclip_5p_count\tsoftclip_3p_count'
              '\tinsertions_5p_count\tinsertions_3p_count\tdeletion_count', file=sys.stderr)
        for contig in self.softclipped_sites:
            for pos in self.softclipped_sites[contig]:
                site = self.softclipped_sites[contig][pos]
                if site.total_count > 0:
                    print('{contig}\t{pos}\t'.format(contig=contig, pos=pos) + str(site), file=sys.stderr)


class SoftclipSite:

    total_count = None

    softclip_5p_reads = None
    softclip_3p_reads = None

    meets_minlength_5p = None
    meets_minlength_3p = None

    insertion_5p_count = None
    insertion_3p_count = None
    runthrough_count = None
    deletion_count = None
    upstream_deletion_count = None
    downstream_deletion_count = None

    keep_softclips_5p = None
    keep_softclips_3p = None

    consensus_sequences_5p = None
    consensus_sequences_3p = None

    def __init__(self):
        self. total_count = 0
        self.softclip_5p_reads = set()
        self.softclip_3p_reads = set()

        self.meets_minlength_5p = False
        self.meets_minlength_3p = False

        self.insertion_5p_count = 0
        self.insertion_3p_count = 0
        self.runthrough_count = 0
        self.deletion_count = 0
        self.upstream_deletion_count = 0
        self.downstream_deletion_count = 0

        self.keep_softclips_5p = True
        self.keep_softclips_3p = True

        self.consensus_sequences_5p = []
        self.consensus_sequences_3p = []

    def add_softclip_5p(self, read, meets_minlength):
        self.softclip_5p_reads.add(read)
        self.total_count += 1

        if not self.meets_minlength_5p and meets_minlength:
            self.meets_minlength_5p = True

    def add_softclip_3p(self, read, meets_minlength):
        self.softclip_3p_reads.add(read)
        self.total_count += 1

        if not self.meets_minlength_3p and meets_minlength:
            self.meets_minlength_3p = True

    def add_runthrough(self):
        self.runthrough_count += 1
        self.total_count += 1

    def add_runthrough_count(self, count):
        self.runthrough_count += count
        self.total_count += count

    def add_insertion_5p(self):
        self.insertion_5p_count += 1
        self.total_count += 1

    def add_insertion_5p_count(self, count):
        self.insertion_5p_count += count
        self.total_count += count

    def add_insertion_3p(self):
        self.insertion_3p_count += 1
        self.total_count += 1

    def add_insertion_3p_count(self, count):
        self.insertion_3p_count += count
        self.total_count += count

    def add_deletion(self):
        self.deletion_count += 1
        self.total_count += 1

    def add_deletion_count(self, count):
        self.deletion_count += count
        self.total_count += count

    def add_upstream_deletion_count(self, count):
        self.upstream_deletion_count += count

    def add_downstream_deletion_count(self, count):
        self.downstream_deletion_count += count

    def get_softclip_5p_count(self):
        return len(self.softclip_5p_reads)

    def get_softclip_3p_count(self):
        return len(self.softclip_3p_reads)

    def get_softclip_ratio_5p(self):
        return len(self.softclip_5p_reads) / self.total_count

    def get_softclip_ratio_3p(self):
        return len(self.softclip_3p_reads) / self.total_count

    def get_upstream_deletion_ratio(self):
        return self.upstream_deletion_count / self.total_count

    def get_downstream_deletion_ratio(self):
        return self.downstream_deletion_count / self.total_count

    def get_softclip_ratio_3p(self):
        return len(self.softclip_3p_reads) / self.total_count

    def get_indel_ratio_5p(self):
        return (self.insertion_5p_count + self.deletion_count) / self.total_count

    def get_indel_ratio_3p(self):
        return (self.insertion_3p_count + self.deletion_count) / self.total_count

    def __str__(self):
        outstring = '{total}\t{runthrough_count}\t{softclip_5p_count}\t{softclip_3p_count}' \
                    '\t{insertions_5p_count}\t{insertions_3p_count}\t{deletion_count}'.format(
            total=self.total_count, runthrough_count=self.runthrough_count,
            softclip_5p_count=len(self.softclip_5p_reads), softclip_3p_count=len(self.softclip_3p_reads),
            insertions_5p_count=self.insertion_5p_count, insertions_3p_count=self.insertion_3p_count,
            deletion_count=self.deletion_count)

        return outstring


class SoftclipConsensus:

    reads = None
    orient = None
    softclip_pos = None
    min_count_consensus = None
    softclipped_seqs = None
    softclipped_qualities = None
    consensus_seqs = None

    def __init__(self, reads, orient, softclip_pos, min_count_consensus):
        self.reads = reads
        self.orient = orient
        self.softclip_pos = softclip_pos
        self.min_count_consensus = min_count_consensus
        self.softclipped_seqs = []
        self.softclipped_qualities = []
        self.consensus_seqs = []

        self.get_softclipped_sequences_and_qualities()
        self.make_consensus()


    def make_consensus(self):

        flanktrie = self.get_flank_sequence_trie(self.softclipped_seqs, self.softclipped_qualities)
        flanktrie_traversal = flanktrie.traverse_all()
        seq_dict = {seq[0]: [seq[1], seq[2]] for seq in flanktrie_traversal}
        seq_clusters = self.get_sequence_clusters(seq_dict)

        consensus_seqs = self.get_cluster_consensus_seqs(seq_clusters)

        if self.orient == '3p':
            consensus_seqs = [seq[::-1] for seq in consensus_seqs]

        merged_counts = self.merge_cluster_counts(seq_clusters, flanktrie)

        self.consensus_seqs = zip(merged_counts, consensus_seqs)


    def get_flank_sequence_trie(self, softclipped_seqs, softclipped_base_qual):

        mytrie = flanktrie.Trie()

        for i in range(len(softclipped_seqs)):
            seq, qual = softclipped_seqs[i], softclipped_base_qual[i]
            mytrie.add(seq.upper(), qual)

        return mytrie


    def get_softclipped_sequences_and_qualities(self):

        self.softclipped_seqs = []
        self.softclipped_qualities = []

        for read in self.reads:

            if self.orient == '3p':
                self.softclipped_seqs.append(sctools.left_softclipped_sequence(read)[::-1])
                self.softclipped_qualities.append(sctools.left_softclip_qualities(read)[::-1])

            elif self.orient == '5p':
                self.softclipped_seqs.append(sctools.right_softclipped_sequence(read))
                self.softclipped_qualities.append(sctools.right_softclip_qualities(read))


    def get_sequence_clusters(self, sequence_dict, perc_similarity=0.75):

        if len(sequence_dict) == 1:
            return [sequence_dict]

        sequences = list(sequence_dict.keys())

        no_match = set(sequences)
        has_match = set()
        pairwise_matches = set()
        for i, j in itertools.combinations(range(len(sequences)), r=2):

            seq1, seq2 = sequences[i], sequences[j]
            minlength = min([len(seq1), len(seq2)])
            dist = levenshtein_distance(seq1[:minlength], seq2[:minlength])
            similarity = 1 - (dist / minlength)

            if similarity >= perc_similarity:
                has_match.add(i)
                has_match.add(j)
                if seq1 in no_match: no_match.remove(seq1)
                if seq2 in no_match: no_match.remove(seq2)
                pairwise_matches.add((i, j))

        pair_order = sorted(list(has_match))
        adjacency_matrix = pd.DataFrame(np.zeros((len(pair_order), len(pair_order))), index=pair_order,
                                        columns=pair_order)

        for pair in pairwise_matches:
            adjacency_matrix.loc[pair[0], pair[1]] = 1
            adjacency_matrix.loc[pair[1], pair[0]] = 1

        n_comp, assignments = connected_components(adjacency_matrix.values)
        final_components = []
        for i in range(n_comp):
            comp = list(np.array(pair_order)[assignments == i])
            final_components.append(comp)

        final_clusters = [{seq: sequence_dict[seq]} for seq in no_match]

        for comp in final_components:
            new_seq_cluster = dict()
            for i in comp:
                new_seq_cluster[sequences[i]] = sequence_dict[sequences[i]]
            final_clusters.append(new_seq_cluster)

        return final_clusters

    def get_cluster_consensus_seqs(self, seq_clusters):
        consensus_seqs = list()

        for clust in seq_clusters:
            seqs = [seq for seq in clust]
            quals = [clust[seq][0] for seq in clust]
            counts = [clust[seq][1] for seq in clust]

            mytrie = flanktrie.Trie()
            mytrie.load_words(seqs, quals, counts)
            consensus = mytrie.make_consensus_word(self.min_count_consensus)

            consensus_seqs.append(consensus)

        return consensus_seqs

    def merge_cluster_counts(self, seq_clusters, flanktrie):

        seq_clusters = [list(cluster.keys()) for cluster in seq_clusters]
        out = []
        for clust in seq_clusters:
            subtrie = flanktrie.make_subtrie(clust)
            total_words = subtrie.total_words
            out.append(total_words)

        return out

    def print_all(self):
        if self.orient == '5p':
            print('SOFTCLIPPED SEQS:')
            for seq in self.softclipped_seqs:
                print('\t'+seq)
            print('CONSENSUS SEQS:')
            for seq in self.consensus_seqs:
                print('\t' + seq[1])
            print()
        else:
            print('SOFTCLIPPED SEQS:')
            for seq in self.softclipped_seqs:
                print('\t' + seq[::-1])
            print('CONSENSUS SEQS:')
            for seq in self.consensus_seqs:
                print('\t' + seq[1])
            print()


if __name__ == '__main__':

    softclip_parser = SoftclipParser(pysam.AlignmentFile(sys.argv[1], 'rb'))
    softclip_parser.parse_softclips()
    softclip_parser.filter_softclips_minlength()
    softclip_parser.filter_softclips_mincount()
    softclip_parser.filter_softclips_mindistance()
    softclip_parser.parse_unclipped_read_info()
    softclip_parser.filter_softclips_count_ratios()
    softclip_parser.filter_softclips_mindistance()
    softclip_parser.make_consensus_sequences()
    softclip_parser.filter_consensus_sequences_minlength()
    softclip_parser.filter_consensus_sequences_mincount()
    softclip_parser.filter_softclips_mindistance()

    print(softclip_parser.make_dataframe().to_string())
