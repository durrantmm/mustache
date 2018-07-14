import sys
from mustache import alignment_tools, output, misc
from mustache.bwa_tools import *
from collections import defaultdict
from Bio import SeqIO
from os.path import join
from snakemake import shell

def get_contigs_with_reads(bam_file):
    contigs = set()
    for r in bam_file.fetch():
        contigs.add(r.reference_name)
    return(list(contigs))

def get_flank_read_count(genome, contig, site, orientation, flank_assembly, unmapped_reads, softclipped_reads,
                         flank_length, outdir, site_name, min_alignment_overlap=1, min_alignment_overlap_count=0):

    inserted_genome = get_inserted_genome(genome, contig, site, orientation, flank_assembly, flank_length)

    tmp_genome = join(outdir, site_name + '.genome.read_count.fasta')
    tmp_fasta = join(outdir, site_name + '.read_count.fasta')
    tmp_alignment = join(outdir, site_name + ".realign.read_count.bam")

    output.write_reads_to_fasta([(site_name + "_ancestral_genome", inserted_genome)], tmp_genome)
    output.write_reads_to_fasta(softclipped_reads + unmapped_reads, tmp_fasta)

    index_genome(tmp_genome, silence=True)
    align_to_genome_single(tmp_fasta, tmp_genome, tmp_alignment)


    total_read_count, junction_read_count = get_inseq_overlapping_read_count(pysam.AlignmentFile(tmp_alignment, 'rb'),
                                                                             len(flank_assembly), orientation,
                                                                             min_alignment_overlap,
                                                                             min_alignment_overlap_count)

    shell("rm {tmp_genome}* {tmp_fasta}* {tmp_alignment}*;")

    return junction_read_count, total_read_count


def get_inseq_overlapping_read_count(bam_file, assembly_length, orientation, min_alignment_overlap, min_alignment_overlap_count):

    unique_all_reads = set()
    unique_junction_reads = set()
    total_overlap_count = 0

    total_reads = 0

    if orientation == 'L':
        for read in bam_file:
            if read.is_unmapped:
                continue

            if read.reference_start <= assembly_length:
                unique_all_reads.add(read.query_name)

            if read.reference_start <= (assembly_length - min_alignment_overlap):
                total_overlap_count += 1

            if read.reference_start <= assembly_length and read.reference_end > assembly_length:
                unique_junction_reads.add(read.query_name)

            total_reads += 1

    elif orientation == 'R':

        contig_length = bam_file.header['SQ'][0]['LN']

        for read in bam_file:
            if read.is_unmapped:
                continue

            if read.reference_end >= (contig_length - assembly_length):
                unique_all_reads.add(read.query_name)

            if read.reference_end >= (contig_length - assembly_length + min_alignment_overlap):
                total_overlap_count += 1

            if read.reference_end >= (contig_length - assembly_length) and \
                            read.reference_start < (contig_length - assembly_length + min_alignment_overlap):
                unique_junction_reads.add(read.query_name)

            total_reads += 1


    total_read_count = len(unique_all_reads)
    junction_read_count = len(unique_junction_reads)

    if total_overlap_count < min_alignment_overlap_count:
        total_read_count = 0
        junction_read_count = 0

    return(total_read_count, junction_read_count)



def get_left_softclipped_reads_at_site(bam_file, contig, left_site, revcomp_left=False):
    left_softclipped_reads = []

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

    return left_softclipped_reads


def get_right_softclipped_reads_at_site(bam_file, contig, right_site, revcomp_left):

    right_softclipped_reads = []
    # First the right softclipped reads.
    for pu in bam_file.pileup(contig, right_site-1, right_site, truncate=True):

        for pr in pu.pileups:
            read = pr.alignment
            if is_right_softclipped(read) and right_softclipped_site(read) == right_site:
                read_tuple = ( read.query_name, read.query_sequence, right_softclipped_sequence(read) )
                right_softclipped_reads.append(read_tuple)

    return right_softclipped_reads


def get_left_unmapped_reads_in_flanks(bam_file, contig, left_site, flank_length, max_mapping_distance):
    left_unmapped_reads = []

    left_flank_start = left_site + 1
    left_flank_end = left_site + 1 + flank_length

    # Now the left site
    for read in bam_file.fetch(contig, left_flank_start, left_flank_end):

        if read.is_reverse and (read.tlen == 0 or abs(read.tlen) > max_mapping_distance):
            read_tuple = (read.query_name, read.get_tag('MT'))
            left_unmapped_reads.append(read_tuple)

    return left_unmapped_reads


def get_right_unmapped_reads_in_flanks(bam_file, contig, right_site, flank_length, max_mapping_distance):

    right_unmapped_reads = []

    right_flank_start = right_site - flank_length
    right_flank_end = right_site

    if right_flank_start < 0:
        right_flank_start = 0

    # First the right site...
    for read in bam_file.fetch(contig, right_flank_start, right_flank_end):

        if not read.is_reverse and (read.tlen == 0 or abs(read.tlen) > max_mapping_distance):
            read_tuple = ( read.query_name, read.get_tag('MT') )
            right_unmapped_reads.append(read_tuple)

    return right_unmapped_reads


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

def is_left_softclipped_strict(read):
    if read.cigartuples[0][0] == 4:
        return True
    else:
        return False

def is_right_softclipped_strict(read):
    if read.cigartuples[-1][0] == 4:
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

def get_right_softclip_length_strict(read):
    if is_right_softclipped_strict(read):
        return read.cigartuples[-1][1]
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

def get_left_softclip_length_strict(read):
    if is_left_softclipped_strict(read):
        return read.cigartuples[0][1]
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

def get_most_common_start_site(bam_file):

    site_counts = defaultdict(int)

    for r in bam_file:
        if r.is_unmapped:
            continue
        if r.reference_start != r.reference_end - 1:
            site_counts[(r.reference_name, r.reference_start)] += 1
            site_counts[(r.reference_name, r.reference_end - 1)] += 1
        else:
            site_counts[(r.reference_name, r.reference_start)] += 1

    try:
        most_common_site = max(list(dict(site_counts).items()), key=lambda x: x[1])[0]
    except ValueError:
        return None

    return most_common_site, site_counts[most_common_site]


def retrieve_flanking_sequence(fasta_file, bam_file, orientation):
    most_common_start_site = get_most_common_start_site(bam_file)

    if not most_common_start_site:
        return None

    start_site = most_common_start_site[0]

    alignment_orientation = get_most_common_orientation_at_site(bam_file, start_site[0], start_site[1])

    fasta = SeqIO.parse(fasta_file, format='fasta')

    match_rec = None
    for rec in fasta:
        if rec.name == start_site[0]:
            match_rec = rec
            break

    if orientation == 'R':
        if alignment_orientation == "FORWARD":
            return str(match_rec.seq)[start_site[1]:].upper()
        elif alignment_orientation == "REVERSE":
            return misc.revcomp(str(match_rec.seq)[:(start_site[1] + 1)]).upper()
        else:
            return None

    elif orientation == "L":
        if alignment_orientation == "FORWARD":
            return misc.revcomp(str(match_rec.seq)[start_site[1]:]).upper()
        elif alignment_orientation == "REVERSE":
            return str(match_rec.seq)[:(start_site[1] + 1)].upper()
        else:
            return None


def retrieve_merged_assembly(assembly, left_bam, right_bam):

    most_common_left_start_site = get_most_common_start_site(left_bam)
    most_common_right_start_site = get_most_common_start_site(right_bam)

    if not most_common_left_start_site or not most_common_right_start_site:
        return None

    left_start_site, left_start_count = most_common_left_start_site
    right_start_site, right_start_count = most_common_right_start_site

    left_contig, left_site = left_start_site
    right_contig, right_site = right_start_site

    if left_contig != right_contig:
        return None
    else:
        fasta = SeqIO.parse(assembly, format='fasta')

        match_contig = None
        for rec in fasta:
            if rec.name == left_contig:
                match_contig = str(rec.seq)
                break

        if left_site < right_site:
            left_orientation = get_most_common_orientation_at_site(left_bam, left_contig, left_site)
            right_orientation = get_most_common_orientation_at_site(right_bam, right_contig, right_site)

            if left_orientation == "FORWARD" and right_orientation == "REVERSE":
                return match_contig[left_site:(right_site+1)]
            else:
                return None

        elif left_site > right_site:
            left_orientation = get_most_common_orientation_at_site(left_bam, left_contig, left_site)
            right_orientation = get_most_common_orientation_at_site(right_bam, right_contig, right_site)

            if left_orientation == "REVERSE" and right_orientation == "FORWARD":
                return match_contig[right_site:(left_site + 1)]
            else:
                return None

        else:
            return None

def get_most_common_orientation_at_site(bam_file, contig, site):
    reverse_count = 0
    forward_count = 0

    total_count = 0
    for pu in bam_file.pileup(contig, site, site + 1, truncate=True):

        for pr in pu.pileups:
            total_count += 1
            read = pr.alignment
            if read.is_reverse:
                reverse_count += 1
            else:
                forward_count += 1

    if forward_count > reverse_count:
        return "FORWARD"
    elif forward_count < reverse_count:
        return "REVERSE"
    else:
        print ("WEIRD READ ORIENTATION ISSUE in bam_tools.py")
        sys.exit()


def get_inserted_genome(genome, contig, site, orientation, flank_assembly, flank_length):

    if orientation == 'L':
        inserted_genome = flank_assembly + genome[(site+1):(site+1+flank_length)]
    elif orientation == 'R':
        inserted_genome = genome[(site-flank_length):site] + flank_assembly

    return inserted_genome


def get_average_region_coverage(bam, contig, start, end):

    running_sum = 0

    for pileupcolumn in bam.pileup(contig, start, end, truncate=True):
        running_sum += pileupcolumn.n

    return running_sum / (end-start)