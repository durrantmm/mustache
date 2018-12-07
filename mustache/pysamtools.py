import warnings
warnings.filterwarnings("ignore")
from mustache.sctools import *

def get_softclipped_seqs_qualities(bam_file, contig, pos, orient, min_alignment_quality, min_alignment_inner_length=21):
    softclip_count = 0
    softclipped_seqs = []
    softclipped_qualities = []

    if orient == 'R':
        start, end = pos - 1, pos
        if start < 0:
            start = 0

    if orient == 'L':
        start, end = pos + 1, pos + 2

        if end > contig_length(bam_file, contig):
            end = contig_length(bam_file, contig)

    for pu in bam_file.pileup(contig, start, end, truncate=True):
        for pr in pu.pileups:
            read = pr.alignment

            if read.mapping_quality < min_alignment_quality:
                continue

            if not read_meets_min_alignment_inner_length(read, min_alignment_inner_length):
                continue

            if orient == 'L':
                if is_left_softclipped_lenient_at_site(read, contig, pos):
                    softclip_count += 1
                    softclipped_seqs.append(left_softclipped_sequence(read))
                    softclipped_qualities.append(left_softclip_qualities(read))

            elif orient == 'R':
                if is_right_softclipped_lenient_at_site(read, contig, pos):
                    softclip_count += 1
                    softclipped_seqs.append(right_softclipped_sequence(read))
                    softclipped_qualities.append(right_softclip_qualities(read))

    if orient == 'L':
        softclipped_seqs = [seq[::-1] for seq in softclipped_seqs]
        softclipped_qualities = [qual[::-1] for qual in softclipped_qualities]

    return softclipped_seqs, softclipped_qualities

def get_left_softclipped_reads_at_site(bam_file, contig, left_site, get_quals=False, softclip_only=False):
    left_softclipped_reads = []
    left_softclipped_quals = []

    start = left_site + 1
    end = left_site + 2

    contig_len = contig_length(bam_file, contig)
    if end > contig_len:
        end = contig_len

    for pu in bam_file.pileup(contig, start, end, truncate=True):

        for pr in pu.pileups:
            read = pr.alignment
            if is_left_softclipped_lenient_at_site(read, contig, left_site):

                if softclip_only:
                    keep_seq = left_softclipped_sequence(read)
                else:
                    keep_seq = read.query_sequence
                    keep_quals = get_query_qualities_ascii(read, bam_file)

                left_softclipped_reads.append(keep_seq)
                left_softclipped_quals.append(keep_quals)

    if get_quals:
        return left_softclipped_reads, left_softclipped_quals
    else:
        return left_softclipped_reads



def get_right_softclipped_reads_at_site(bam_file, contig, right_site, get_quals=False, softclip_only=False):
    right_softclipped_reads = []
    right_softclipped_quals = []

    start = right_site - 1
    end = right_site

    if start < 0:
        start = 0

    for pu in bam_file.pileup(contig, start, end, truncate=True):

        for pr in pu.pileups:
            read = pr.alignment
            if is_right_softclipped_lenient_at_site(read, contig, right_site):
                if softclip_only:
                    keep_seq = right_softclipped_sequence(read)
                else:
                    keep_seq = read.query_sequence
                    keep_quals = get_query_qualities_ascii(read, bam_file)

                right_softclipped_reads.append(keep_seq)
                right_softclipped_quals.append(keep_quals)

    if get_quals:
        return right_softclipped_reads, right_softclipped_quals
    else:
        right_softclipped_reads


def get_right_unmapped_reads(bam_file, contig, right_site, get_quals=False, search_region_length=500):
    right_unmapped_reads = []
    right_unmapped_quals = []

    start = right_site - search_region_length
    end = right_site

    if start < 0:
        start = 0

    for read in bam_file.fetch(contig, start, end):

        if (not read.is_reverse) and read.mate_is_unmapped:
            right_unmapped_reads.append(read.get_tag('MT'))
            right_unmapped_quals.append(read.get_tag('MQ'))

    if get_quals:
        return right_unmapped_reads, right_unmapped_quals
    else:
        return right_unmapped_reads


def get_left_unmapped_reads(bam_file, contig, left_site, get_quals=False, search_region_length=500):
    left_unmapped_reads = []
    left_unmapped_quals = []

    start = left_site + 1
    end = left_site + 1 + search_region_length

    contig_len = contig_length(bam_file, contig)
    if end > contig_len:
        end = contig_len

    for read in bam_file.fetch(contig, start, end):

        if read.is_reverse and read.mate_is_unmapped:
            left_unmapped_reads.append(read.get_tag('MT'))
            left_unmapped_quals.append(read.get_tag('MQ'))

    if get_quals:
        return left_unmapped_reads, left_unmapped_quals
    else:
        return left_unmapped_reads


def contig_length(bam, contig):
    return dict(zip(bam.references, bam.lengths))[contig]

def count_runthrough_reads(bam, contig, site, min_qual=20, min_alignment_inner_length=21):

    if site  < 0 or site >= contig_length(bam, contig):
        return 0
    else:
        count = 0
        for read in bam.fetch(contig, site, site+1):

            if read.mapping_quality < min_qual:
                continue

            if not read_meets_min_alignment_inner_length(read, min_alignment_inner_length):
                continue

            count += 1

        return count

def count_softclipped_reads(bam, contig, site, min_qual=20, min_alignment_inner_length=21):

    if site  < 0 or site >= contig_length(bam, contig):
        return 0
    else:
        count = 0
        for read in bam.fetch(contig, site+1, site+2):

            if read.mapping_quality < min_qual:
                continue

            if not read_meets_min_alignment_inner_length(read, min_alignment_inner_length):
                continue

            if is_left_softclipped_lenient_at_site(read, contig, site):
                count += 1

        for read in bam.fetch(contig, site - 1, site):

            if read.mapping_quality < min_qual:
                continue

            if not read_meets_min_alignment_inner_length(read, min_alignment_inner_length):
                continue

            if is_right_softclipped_lenient_at_site(read, contig, site):
                count += 1


        return count

def get_query_qualities_ascii(read, bam):
    return read.tostring(bam).split('\t')[10]

def query_qualities_to_phred(quals):
    return [ord(q)-33 for q in quals]

def get_perc_identity(read):
    total = read.query_alignment_length
    matches = len(read.get_aligned_pairs(matches_only=True))
    identity = matches/total
    return(identity)
