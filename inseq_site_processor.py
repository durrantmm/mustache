import sys
import pysam
import flank_assembler, alignment_tools
import misc
import output
from collections import OrderedDict
import click
import pandas as pd
pd.set_option('display.height', 1000)
pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)


def process_candidate_site(bam_file, contig, site, flank_length, max_mapping_distance, outdir, output_prefix):
    out_dict = OrderedDict([
        ('contig', [contig]),
        ('left_site', [site[0]]),
        ('right_site', [site[1]]),
        ('flank_length', [flank_length]),
        ('direct_inseq_read_count', None),
        ('direct_ancestral_read_count', None),
        ('direct_inseq_freq', None),
        ('flanking_inseq_read_count', None),
        ('flanking_ancestral_read_count', None),
        ('flanking_inseq_freq', None),
        ('right_assembly_length', None),
        ('left_assembly_length', None),
        ('right_assembly', None),
        ('left_assembly', None)
    ])

    site_name = contig + ':' + str(site[0]) + '-' + str(site[1])

    click.echo('Processing the site %s' % site_name)
    bam_path = bam_file.filename.decode('UTF-8')
    bam_file = pysam.AlignmentFile(bam_path, 'rb')



    # Beginning assembly
    click.echo('\tPerforming assembly...')
    left_softclipped_reads, right_softclipped_reads = get_softclipped_reads_at_site(bam_file, contig, site, revcomp_left=True)
    left_unmapped_reads, right_unmapped_reads = unmapped_reads_in_flanks(bam_file, contig, site, flank_length, max_mapping_distance)

    # Assembling the right flank
    right_inseq_assembly = flank_assembler.assemble_flank(right_softclipped_reads, right_unmapped_reads)
    out_dict['right_assembly'] = [right_inseq_assembly]
    out_dict['right_assembly_length'] = [len(right_inseq_assembly)]
    click.echo('\t\tRight assembly:')
    click.echo(misc.wrap_string(right_inseq_assembly, newline_char='\t\t'))
    click.echo()

    # Assembling the left flank
    left_inseq_assembly = misc.revcomp(flank_assembler.assemble_flank(left_softclipped_reads, left_unmapped_reads))
    out_dict['left_assembly'] = [left_inseq_assembly]
    out_dict['left_assembly_length'] = [len(left_inseq_assembly)]
    click.echo('\t\tLeft assembly:')
    click.echo(misc.wrap_string(left_inseq_assembly, newline_char='\t\t'))
    click.echo()

    # Now get direct insertion sequence reads
    click.echo('\tGetting IS overlapping reads directly at the insertion sequence site...')
    direct_inseq_reads = get_direct_inseq_reads(bam_file, contig, site, right_inseq_assembly, left_inseq_assembly)
    out_dict['direct_inseq_read_count'] = [len(direct_inseq_reads)]
    click.echo('\t\t%d IS overlapping reads directly at insertion site' % len(direct_inseq_reads))
    click.echo('\t\tWriting direct IS reads to file...')
    output.write_site_reads(direct_inseq_reads, bam_file, outdir, output_prefix, site_name+'.direct_inseq')
    del direct_inseq_reads
    click.echo()

    # Now get unique site ancestral reads
    click.echo('\tGetting ancestral reads directly overlapping site...')
    direct_ancestral_reads = get_direct_ancestral_reads(bam_file, contig, site, right_inseq_assembly, left_inseq_assembly)
    out_dict['direct_ancestral_read_count'] = [len(direct_ancestral_reads)]
    click.echo('\t\t%d ancestral reads directly at insertion site' % len(direct_ancestral_reads))
    click.echo('\t\tWriting direct ancestral reads to file...')
    output.write_site_reads(direct_ancestral_reads, bam_file, outdir, output_prefix, site_name + '.direct_ancestral')
    del direct_ancestral_reads
    click.echo()

    # Quickly calculate the direct inseq frequency
    click.echo('\tCalculating the direct IS frequency...')
    total_direct_reads = float(out_dict['direct_inseq_read_count'][0] + out_dict['direct_ancestral_read_count'][0])
    direct_inseq_freq = out_dict['direct_inseq_read_count'][0] / total_direct_reads
    out_dict['direct_inseq_freq'] = [direct_inseq_freq]
    click.echo('\t\tThe direct IS frequency is %f' % direct_inseq_freq)
    click.echo()

    # Now get reads that flank the insertion site and map to the insertion
    click.echo('\tGetting reads that flank the IS site, and map to the relevant sequence...')
    flanking_inseq_reads = get_flanking_inseq_reads(bam_file, contig, site, flank_length, right_inseq_assembly, left_inseq_assembly, max_mapping_distance)
    out_dict['flanking_inseq_read_count'] = [len(flanking_inseq_reads)]
    click.echo('\t\t%d reads flanking the insertion site' % len(flanking_inseq_reads))
    click.echo('\t\tWriting flanking IS reads to file...')
    output.write_site_reads(flanking_inseq_reads, bam_file, outdir, output_prefix, site_name + '.flanking_inseq')
    del flanking_inseq_reads
    click.echo()

    # Now get the ancestral reads that flank the insertion...
    click.echo('\tGetting ancestral reads that flank the IS site...')
    flanking_ancestral_reads = get_flanking_ancestral_reads(bam_file, contig, site, flank_length, right_inseq_assembly,
                                                            left_inseq_assembly, max_mapping_distance)
    out_dict['flanking_ancestral_read_count'] = [len(flanking_ancestral_reads)]
    click.echo('\t\t%d ancestral reads flanking the insertion site' % len(flanking_ancestral_reads))
    click.echo('\t\tWriting flanking ancestral reads to file...')
    output.write_site_reads(flanking_ancestral_reads, bam_file, outdir, output_prefix, site_name + '.flanking_ancestral')
    del flanking_ancestral_reads
    click.echo()

    # Quickly calculate the flanking inseq frequency
    click.echo('\tCalculating the flanking IS frequency...')
    total_direct_reads = float(out_dict['flanking_inseq_read_count'][0] + out_dict['flanking_ancestral_read_count'][0])
    flanking_inseq_freq = out_dict['flanking_inseq_read_count'][0] / total_direct_reads
    out_dict['flanking_inseq_freq'] = [flanking_inseq_freq]
    click.echo('\t\tThe flanking IS frequency is %f' % flanking_inseq_freq)
    click.echo()

    return pd.DataFrame(out_dict)


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

    # First the right site...
    for read in bam_file.fetch(contig, right_site-flank_length, right_site):

        if not read.is_reverse and (read.tlen == 0 or abs(read.tlen) > max_mapping_distance):
            read_tuple = ( read.query_name, read.get_tag('MT') )
            right_unmapped_reads.append(read_tuple)

    # Now the left site
    for read in bam_file.fetch(contig, left_site+1, left_site+1+flank_length):

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

    left_flank_start = left_site + 1
    left_flank_end = left_site + 1 + flank_length

    right_assembly_mapper = alignment_tools.QuickMapper(right_assembly)
    left_assembly_mapper = alignment_tools.QuickMapper(left_assembly)

    # First the right site...
    # Starting with softclipped reads...
    for pu in bam_file.pileup(contig, right_site - 1, right_site, truncate=True):

        for pr in pu.pileups:
            read = pr.alignment

            if is_right_softclipped(read):

                site_clipped_read = clip_read_at_right_site(read, right_site)
                align_score = alignment_tools.left_alignment_score(site_clipped_read, right_assembly)

                if align_score > 0:

                    if read.is_reverse and not (read.tlen == 0 or abs(read.tlen) > max_mapping_distance) and \
                        right_flank_start <= read.next_reference_start < right_flank_end and \
                        read.query_name not in keep_read_names:

                        read_mate = bam_file.mate(read)
                        keep_read_names.add(read_mate.query_name)
                        flanking_reads.append(read_mate)

                    elif (not read.is_reverse) and (read.tlen == 0 or abs(read.tlen) > max_mapping_distance) and \
                        right_assembly_mapper.maps_to_reference(read.get_tag('MT')) and \
                        read.query_name not in keep_read_names:

                        keep_read_names.add(read.query_name)
                        flanking_reads.append(read)

    # Now the unmapped reads at the right site...
    for read in bam_file.fetch(contig, right_flank_start, right_flank_end):

        if not read.is_reverse and (read.tlen == 0 or abs(read.tlen) > max_mapping_distance) and \
            right_assembly_mapper.maps_to_reference(read.get_tag('MT')) and \
            read.query_name not in keep_read_names:

            keep_read_names.add(read.query_name)
            flanking_reads.append(read)

    # Now the left site
    # Starting with softclipped reads
    for pu in bam_file.pileup(contig, left_site + 1, left_site + 2, truncate=True):

        for pr in pu.pileups:
            read = pr.alignment

            if is_left_softclipped(read):

                site_clipped_read = clip_read_at_left_site(read, left_site)
                align_score = alignment_tools.right_alignment_score(site_clipped_read, left_assembly)

                if align_score > 0 and read.query_name not in keep_read_names:

                    if not read.is_reverse and not (read.tlen == 0 or abs(read.tlen) > max_mapping_distance) and \
                        left_flank_start <= read.next_reference_start < left_flank_end and \
                        read.query_name not in keep_read_names:

                        read_mate = bam_file.mate(read)
                        keep_read_names.add(read_mate.query_name)
                        flanking_reads.append(read_mate)

                    elif read.is_reverse and (read.tlen == 0 or abs(read.tlen) > max_mapping_distance) and \
                        left_assembly_mapper.maps_to_reference(read.get_tag('MT')) and \
                        read.query_name not in keep_read_names:

                        keep_read_names.add(read.query_name)
                        flanking_reads.append(read)

    # Now the unmapped reads at the left site...
    for read in bam_file.fetch(contig, left_flank_start, left_flank_end):

        if read.is_reverse and (read.tlen == 0 or abs(read.tlen) > max_mapping_distance) and \
            left_assembly_mapper.maps_to_reference(read.get_tag('MT')) and \
            read.query_name not in keep_read_names:

            keep_read_names.add(read.query_name)
            flanking_reads.append(read)

    return flanking_reads

def get_flanking_ancestral_reads(bam_file, contig, site, flank_length, right_assembly, left_assembly, max_mapping_distance):
    keep_read_names = set()
    flanking_reads = []

    left_site = site[0]
    right_site = site[1]

    right_flank_start = right_site - flank_length
    right_flank_end = right_site

    left_flank_start = left_site + 1
    left_flank_end = left_site + 1 + flank_length

    # First the reads oriented toward the right site...
    for read in bam_file.fetch(contig, right_flank_start, right_flank_end):

        if (not read.is_reverse) and read.tlen != 0 and abs(read.tlen) <= max_mapping_distance and \
            read.reference_start + read.tlen >= right_flank_end:

            read_mate = bam_file.mate(read)

            if is_right_softclipped(read):

                site_clipped_read = clip_read_at_right_site(read, right_site)
                align_score = alignment_tools.left_alignment_score(site_clipped_read, right_assembly)

                if align_score > 0:
                    continue

            if is_right_softclipped(read_mate):

                site_clipped_read = clip_read_at_right_site(read_mate, right_site)
                align_score = alignment_tools.left_alignment_score(site_clipped_read, right_assembly)

                if align_score > 0:
                    continue

            if is_left_softclipped(read):

                site_clipped_read = clip_read_at_left_site(read, left_site)
                align_score = alignment_tools.right_alignment_score(site_clipped_read, left_assembly)

                if align_score > 0:
                    continue

            if is_left_softclipped(read_mate):

                site_clipped_read = clip_read_at_left_site(read_mate, left_site)
                align_score = alignment_tools.right_alignment_score(site_clipped_read, left_assembly)

                if align_score > 0:
                    continue

            if read.query_name not in keep_read_names:
                keep_read_names.add(read.query_name)
                flanking_reads.append(read)
                flanking_reads.append(read_mate)

    # Now the reads at the left site...
    for read in bam_file.fetch(contig, left_flank_start, left_flank_end):

        if read.is_reverse and read.tlen != 0 and abs(read.tlen) <= max_mapping_distance and \
            read.reference_end + read.tlen <= left_flank_start:

            read_mate = bam_file.mate(read)

            if is_right_softclipped(read):

                site_clipped_read = clip_read_at_right_site(read, right_site)
                align_score = alignment_tools.left_alignment_score(site_clipped_read, right_assembly)

                if align_score > 0:
                    continue

            if is_right_softclipped(read_mate):

                site_clipped_read = clip_read_at_right_site(read_mate, right_site)
                align_score = alignment_tools.left_alignment_score(site_clipped_read, right_assembly)

                if align_score > 0:
                    continue

            if is_left_softclipped(read):

                site_clipped_read = clip_read_at_left_site(read, left_site)
                align_score = alignment_tools.right_alignment_score(site_clipped_read, left_assembly)

                if align_score > 0:
                    continue

            if is_left_softclipped(read_mate):

                site_clipped_read = clip_read_at_left_site(read_mate, left_site)
                align_score = alignment_tools.right_alignment_score(site_clipped_read, left_assembly)

                if align_score > 0:
                    continue

            if read.query_name not in keep_read_names:
                keep_read_names.add(read.query_name)
                flanking_reads.append(read)
                flanking_reads.append(read_mate)

    return flanking_reads


def is_right_softclipped(read):
    return read.cigartuples[-1][0] == 4

def is_left_softclipped(read):
    return read.cigartuples[0][0] == 4

def right_softclipped_site(read):
    return read.get_reference_positions()[-1] + 1

def left_softclipped_site(read):
    return read.get_reference_positions()[0] - 1

def right_softclipped_sequence(read):
    if is_right_softclipped(read):
        return read.query_sequence[-read.cigartuples[-1][1]:]
    else:
        return ''

def left_softclipped_sequence(read):
    return read.query_sequence[:read.cigartuples[0][1]]

def clip_read_at_right_site(read, right_site):
    clip_length = read.cigartuples[-1][1] + read.get_reference_positions()[-1] - right_site + 1
    clipped_read = read.query_sequence[-clip_length:]
    return clipped_read

def clip_read_at_left_site(read, left_site):
    clip_length = read.cigartuples[0][1] + left_site - read.get_reference_positions()[0] + 1
    clipped_read = read.query_sequence[:clip_length]
    return clipped_read

def read_runs_through_site(read, left_site, right_site):
    return read.get_reference_positions()[0] <= left_site and read.get_reference_positions()[-1] >= right_site