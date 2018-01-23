import sys
import alignment_tools
import misc

def get_softclipped_reads_at_site(bam_file, contig, site, revcomp_left=False):

    left_softclipped_reads = []
    right_softclipped_reads = []

    left_site = site[0]
    right_site = site[1]

    # First the right softclipped reads.
    for pu in bam_file.pileup(contig, right_site-1, right_site, truncate=True):

        for pr in pu.pileups:
            read = pr.alignment
            if is_right_softclipped(read) and right_softclipped_site(read) == right_site:
                read_tuple = ( read.query_name, read.query_sequence, right_softclipped_sequence(read) )
                right_softclipped_reads.append(read_tuple)

    # Now the left softclipped reads:
    softclipped_reads = []

    left_site = site[0]
    for pu in bam_file.pileup(contig, left_site + 1, left_site + 2, truncate=True):

        for pr in pu.pileups:
            read = pr.alignment
            if is_left_softclipped(read) and left_softclipped_site(read) == left_site:
                if not revcomp_left:
                    read_softclip_seq = left_softclipped_sequence(read)
                    read_full_seq = read.query_sequence
                else:
                    read_softclip_seq = misc.revcomp(left_softclipped_sequence(read))
                    read_full_seq = misc.revcomp(read.query_sequence)

                read_tuple = (read.query_name, read_full_seq, read_softclip_seq)
                left_softclipped_reads.append(read_tuple)

    return left_softclipped_reads, right_softclipped_reads


def unmapped_reads_in_flanks(bam_file, contig, site, flank_length, max_mapping_distance):
    left_unmapped_reads = []
    right_unmapped_reads = []

    left_site = site[0]
    right_site = site[1]

    right_flank_start = right_site - flank_length
    right_flank_end = right_site

    if right_flank_start < 0:
        right_flank_start = 0

    left_flank_start = left_site + 1
    left_flank_end = left_site + 1 + flank_length

    # First the right site...
    for read in bam_file.fetch(contig, right_flank_start, right_flank_end):

        if not read.is_reverse and (read.tlen == 0 or abs(read.tlen) > max_mapping_distance):
            read_tuple = ( read.query_name, read.get_tag('MT') )
            right_unmapped_reads.append(read_tuple)

    # Now the left site
    for read in bam_file.fetch(contig, left_flank_start, left_flank_end):

        if read.is_reverse and (read.tlen == 0 or abs(read.tlen) > max_mapping_distance):
            read_tuple = ( read.query_name, read.get_tag('MT') )
            left_unmapped_reads.append(read_tuple)


    return left_unmapped_reads, right_unmapped_reads


def get_direct_inseq_reads(bam_file, contig, site, right_assembly, left_assembly):

    keep_read_names = set()
    softclipped_reads = []

    left_site = site[0]
    right_site = site[1]

    # First the right softclipped reads
    for pu in bam_file.pileup(contig, right_site-1, right_site, truncate=True):
        for pr in pu.pileups:
            read = pr.alignment

            if is_right_softclipped(read):

                site_clipped_read = clip_read_at_right_site(read, right_site)
                align_score = alignment_tools.left_alignment_score(site_clipped_read, right_assembly)

                if align_score > 0 and read.query_name not in keep_read_names:
                    softclipped_reads.append(read)
                    keep_read_names.add(read.query_name)

    # Now the left softclipped reads.
    for pu in bam_file.pileup(contig, left_site + 1, left_site + 2, truncate=True):
        for pr in pu.pileups:
            read = pr.alignment

            if is_left_softclipped(read):

                site_clipped_read = clip_read_at_left_site(read, left_site)
                align_score = alignment_tools.right_alignment_score(site_clipped_read, left_assembly)

                if align_score > 0 and read.query_name not in keep_read_names:
                    softclipped_reads.append(read)
                    keep_read_names.add(read.query_name)

    return softclipped_reads


def get_direct_ancestral_reads(bam_file, contig, site, right_assembly, left_assembly):

    keep_read_names = set()
    ancestral_reads = []

    # First we'll try rightmost non-inseq sight
    left_site = site[0]
    right_site = site[1]

    for pu in bam_file.pileup(contig, right_site-1, right_site, truncate=True):
        for pr in pu.pileups:
            read = pr.alignment

            if read_runs_through_site(read, left_site, right_site):

                if is_right_softclipped(read):

                    site_clipped_read = clip_read_at_right_site(read, right_site)
                    align_score = alignment_tools.left_alignment_score(site_clipped_read, right_assembly)
                    if align_score > 0 and read.query_name not in keep_read_names:
                        continue

                if is_left_softclipped(read):

                    site_clipped_read = clip_read_at_left_site(read, left_site)
                    align_score = alignment_tools.right_alignment_score(site_clipped_read, left_assembly)

                    if align_score > 0 and read.query_name not in keep_read_names:
                        continue

                if read.query_name not in keep_read_names:
                    keep_read_names.add(read.query_name)
                    ancestral_reads.append(read)

    # Now leftmost non-inseq site
    left_site = site[0]
    right_site = site[1]
    for pu in bam_file.pileup(contig, left_site + 1, left_site+2, truncate=True):
        for pr in pu.pileups:
            read = pr.alignment

            if read_runs_through_site(read, left_site, right_site):

                if is_right_softclipped(read):

                    site_clipped_read = clip_read_at_right_site(read, right_site)
                    align_score = alignment_tools.left_alignment_score(site_clipped_read, right_assembly)

                    if align_score > 0 and read.query_name not in keep_read_names:
                        continue

                if is_left_softclipped(read):

                    site_clipped_read = clip_read_at_left_site(read, left_site)
                    align_score = alignment_tools.right_alignment_score(site_clipped_read, left_assembly)

                    if align_score > 0 and read.query_name not in keep_read_names:
                        continue

                if read.query_name not in keep_read_names:
                    keep_read_names.add(read.query_name)
                    ancestral_reads.append(read)

    return ancestral_reads

def get_flanking_inseq_reads(bam_file, contig, site, flank_length, right_assembly, left_assembly, max_mapping_distance):

    keep_read_names = set()
    flanking_reads = []

    left_site = site[0]
    right_site = site[1]

    right_flank_start = right_site - flank_length
    right_flank_end = right_site

    if right_flank_start < 0:
        right_flank_start = 0

    left_flank_start = left_site + 1
    left_flank_end = left_site + 1 + flank_length

    right_assembly_mapper = alignment_tools.QuickMapper(right_assembly)
    left_assembly_mapper = alignment_tools.QuickMapper(left_assembly)

    # First the right site...
    for read in bam_file.fetch(contig, right_flank_start, right_flank_end):

        if not read.is_reverse and (read.tlen == 0 or abs(read.tlen) > max_mapping_distance) and \
            right_assembly_mapper.maps_within_reference(read.get_tag('MT')) and \
            read.query_name not in keep_read_names and read.reference_end-1 <= left_site :

            if is_right_softclipped(read):

                site_clipped_read = clip_read_at_right_site(read, right_site)
                align_score = alignment_tools.left_alignment_score(site_clipped_read, right_assembly)

                if align_score > 0:
                    # Excluding softclipped reads
                    continue
                else:
                    keep_read_names.add(read.query_name)
                    flanking_reads.append(read)
            else:

                keep_read_names.add(read.query_name)
                flanking_reads.append(read)

       # Now the left site
    for read in bam_file.fetch(contig, left_flank_start, left_flank_end):

        if read.is_reverse and (read.tlen == 0 or abs(read.tlen) > max_mapping_distance) and \
            left_assembly_mapper.maps_within_reference(read.get_tag('MT')) and \
            read.query_name not in keep_read_names and read.reference_start >= right_site:

            if is_left_softclipped(read):

                site_clipped_read = clip_read_at_left_site(read, left_site)
                align_score = alignment_tools.right_alignment_score(site_clipped_read, left_assembly)

                if align_score > 0:
                    # Excluding softclipped_reads
                    continue
                else:
                    keep_read_names.add(read.query_name)
                    flanking_reads.append(read)
            else:
                keep_read_names.add(read.query_name)
                flanking_reads.append(read)

    return flanking_reads


def get_flanking_ancestral_reads(bam_file, contig, site, flank_length, right_assembly, left_assembly, max_mapping_distance):

    forward_keep_read_names = set()
    reverse_keep_read_names = set()

    flanking_reads = []

    left_site = site[0]
    right_site = site[1]

    right_flank_start = right_site - flank_length
    right_flank_end = right_site

    if right_flank_start < 0:
        right_flank_start = 0

    left_flank_start = left_site + 1
    left_flank_end = left_site + 1 + flank_length

    # First the reads oriented toward the right site...
    for read in bam_file.fetch(contig, right_flank_start, right_flank_end):

        if (not read.is_reverse) and read.tlen != 0 and abs(read.tlen) <= max_mapping_distance and \
            read.next_reference_start >= right_site and read.reference_end-1 <= left_site:

            if is_right_softclipped(read):

                site_clipped_read = clip_read_at_right_site(read, right_site)
                align_score = alignment_tools.left_alignment_score(site_clipped_read, right_assembly)

                if align_score > 0:
                    continue

            if is_left_softclipped(read):

                site_clipped_read = clip_read_at_left_site(read, left_site)
                align_score = alignment_tools.right_alignment_score(site_clipped_read, left_assembly)

                if align_score > 0:
                    continue

            if read.query_name not in forward_keep_read_names:
                forward_keep_read_names.add(read.query_name)
                flanking_reads.append(read)

    for read in bam_file.fetch(contig, left_flank_start, left_flank_end):

        if read.is_reverse and read.query_name in forward_keep_read_names:
            reverse_keep_read_names.add(read.query_name)
            flanking_reads.append(read)

    final_flanking_reads = [read for read in flanking_reads if read.query_name in forward_keep_read_names & reverse_keep_read_names]

    return final_flanking_reads

def calculate_direct_inseq_scaling_factor(n_inseq_reads, average_read_length, max_softclip_length):
    numerator = float(n_inseq_reads)
    denominator = float(4*average_read_length - 4*(average_read_length-max_softclip_length))
    scaling_factor =  numerator / denominator
    return scaling_factor


def calculate_direct_ancestral_scaling_factor(n_inseq_reads, average_read_length, direct_repeat_length):
    numerator = float(n_inseq_reads)
    denominator = 2*(average_read_length - direct_repeat_length - 1)
    scaling_factor =  numerator / denominator

    return scaling_factor

def is_right_softclipped(read):
    if read.cigartuples[-1][0] == 4:
        return True
    elif read.query_sequence[-1] != read.get_reference_sequence()[-1]:
        return True
    else:
        return False

def is_left_softclipped(read):
    if read.cigartuples[0][0] == 4:
        return True
    elif read.query_sequence[0] != read.get_reference_sequence()[0]:
        return True
    else:
        return False

def get_right_softclip_length(read):
    if is_right_softclipped(read):
        if read.cigartuples[-1][0] == 4:
            return read.cigartuples[-1][1]
        else:
            return 1
    else:
        return 0

def get_left_softclip_length(read):
    if is_left_softclipped(read):
        if read.cigartuples[0][0] == 4:
            return read.cigartuples[0][1]
        else:
            return 1
    else:
        return 0

def right_softclipped_site(read):
    if read.cigartuples[-1][0] == 4:
        return read.get_reference_positions()[-1] + 1
    elif read.query_sequence[-1] != read.get_reference_sequence()[-1]:
        return read.get_reference_positions()[-1]

def left_softclipped_site(read):
    if read.cigartuples[0][0] == 4:
        return read.get_reference_positions()[0] - 1
    elif read.query_sequence[0] != read.get_reference_sequence()[0]:
        return read.get_reference_positions()[0]

def right_softclipped_sequence(read):
    if is_right_softclipped(read):
        if read.cigartuples[-1][0] == 4:
            return read.query_sequence[-read.cigartuples[-1][1]:]
        else:
            return read.query_sequence[-1]
    else:
        return ''

def left_softclipped_sequence(read):
    if is_left_softclipped(read):
        if read.cigartuples[0][0] == 4:
            return read.query_sequence[:read.cigartuples[0][1]]
        else:
            return read.query_sequence[0]
    else:
        return ''

def get_last_unclipped_reference_position(read):
    if is_right_softclipped(read):
        if read.cigartuples[-1][0] == 4:
            return read.reference_end - 1
        else:
            return read.reference_end - 2
    else:
        return read.reference_end - 1

def get_first_unclipped_reference_position(read):
    if is_left_softclipped(read):
        if read.cigartuples[0][0] == 4:
            return read.reference_start
        else:
            return read.reference_start + 1
    else:
        return read.reference_start

def clip_read_at_right_site(read, right_site):
    clip_length = get_right_softclip_length(read) + get_last_unclipped_reference_position(read) - right_site + 1
    clipped_read = read.query_sequence[-clip_length:]
    return clipped_read

def clip_read_at_left_site(read, left_site):
    clip_length = get_left_softclip_length(read) + left_site - get_first_unclipped_reference_position(read) + 1
    clipped_read = read.query_sequence[:clip_length]
    return clipped_read

def read_runs_through_site(read, left_site, right_site):
    return get_first_unclipped_reference_position(read) <= left_site and get_last_unclipped_reference_position(read) >= right_site