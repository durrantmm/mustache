import sys
import pysam
import flank_assembler, alignment_tools
import misc
import output
from collections import OrderedDict
import click
import numpy as np
import pandas as pd
from os.path import join
from bam_ops import *

pd.set_option('display.max_rows', 1000)
pd.set_option('display.max_columns', 1000)
pd.set_option('display.width', 1000)\


def process_candidate_site(bam_file, contig, site, flank_length, max_mapping_distance, max_softclip_length, average_read_length, outdir, output_prefix):
    out_dict = OrderedDict([
        ('contig', [contig]),
        ('left_site', [site[0]]),
        ('right_site', [site[1]]),
        ('direct_repeat_length', [site[1]-site[0]-1]),
        ('flank_length', [flank_length]),
        ('direct_inseq_read_count', np.nan),
        ('direct_ancestral_read_count', np.nan),
        ('direct_inseq_freq', np.nan),
        ('adj_direct_inseq_freq', np.nan),
        ('flanking_inseq_read_count', np.nan),
        ('flanking_ancestral_read_count', np.nan),
        ('flanking_inseq_freq', np.nan),
        ('right_assembly_length', np.nan),
        ('left_assembly_length', np.nan),
        ('merged_assembly_length', np.nan),
        ('right_assembly', np.nan),
        ('left_assembly', np.nan),
        ('merged_assembly', np.nan)
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

    # Now merge the assemblies if possible...
    right_inseq_assembly, left_inseq_assembly, merged_assembly = alignment_tools.merge_flank_assemblies(right_inseq_assembly, left_inseq_assembly)
    if merged_assembly:
        click.echo('\t\tMerged assembly:')
        click.echo(misc.wrap_string(merged_assembly, newline_char='\t\t'))
        out_dict['merged_assembly'] = [merged_assembly]
        out_dict['merged_assembly_length'] = [len(merged_assembly)]
        right_inseq_assembly = merged_assembly
        left_inseq_assembly = merged_assembly
    else:
        click.echo('\t\tThe assembly did not merge...')
    click.echo()

    # Now get direct insertion sequence reads
    direct_inseq_reads = hndl_get_direct_inseq_reads(bam_file, contig, site, right_inseq_assembly, left_inseq_assembly,
                                                     output_prefix, outdir, site_name, out_dict)

    # Now get unique site ancestral reads
    direct_ancestral_reads = hndl_get_ancestral_inseq_reads(bam_file, contig, site, right_inseq_assembly,
                                                            left_inseq_assembly, output_prefix,outdir, site_name,
                                                            out_dict)

    # Estimate the direct inseq frequency
    hndl_estimate_direct_inseq_freq(len(direct_inseq_reads), len(direct_ancestral_reads), out_dict)

    # Estimate the adjusted direct inseq frequency
    hndl_estimate_adj_direct_inseq_freq(max_softclip_length, average_read_length, len(direct_inseq_reads),
                                        len(direct_ancestral_reads), out_dict)

    del direct_ancestral_reads
    del direct_inseq_reads

    # Now get reads that flank the insertion site and map to the insertion
    flanking_inseq_reads = hndl_get_flanking_inseq_reads(bam_file, contig, site, flank_length,
                                                         right_inseq_assembly, left_inseq_assembly,
                                                         max_mapping_distance, output_prefix, outdir,
                                                         site_name, out_dict)

    # Now get the ancestral reads that flank the insertion...
    flanking_ancestral_reads = hndl_get_flanking_ancestral_reads(bam_file, contig, site, flank_length, right_inseq_assembly,
                                                                 left_inseq_assembly, max_mapping_distance,
                                                                 output_prefix, outdir, site_name, out_dict)

    # Quickly calculate the flanking inseq frequency
    hndl_estimate_flanking_inseq_freq(len(flanking_inseq_reads), len(flanking_ancestral_reads), out_dict)

    del flanking_inseq_reads
    del flanking_ancestral_reads

    return pd.DataFrame(out_dict)


def process_call_site(bam_file, contig, site, right_assembly, left_assembly, merged_assembly, flank_length,
                      max_mapping_distance, max_softclip_length, average_read_length, outdir, output_prefix):

    out_dict = OrderedDict([
        ('contig', [contig]),
        ('left_site', [site[0]]),
        ('right_site', [site[1]]),
        ('direct_repeat_length', [site[1] - site[0] - 1]),
        ('flank_length', [flank_length]),
        ('direct_inseq_read_count', np.nan),
        ('direct_ancestral_read_count', np.nan),
        ('direct_inseq_freq', np.nan),
        ('adj_direct_inseq_freq', np.nan),
        ('flanking_inseq_read_count', np.nan),
        ('flanking_ancestral_read_count', np.nan),
        ('flanking_inseq_freq', np.nan),
        ('right_assembly_length', len(right_assembly)),
        ('left_assembly_length', len(left_assembly)),
        ('merged_assembly_length', (lambda: len(merged_assembly) if isinstance(merged_assembly, str) else np.nan)()),
        ('right_assembly', right_assembly),
        ('left_assembly', left_assembly),
        ('merged_assembly', merged_assembly)
    ])

    site_name = contig + ':' + str(site[0]) + '-' + str(site[1])

    click.echo('Processing the site %s' % site_name)
    bam_path = bam_file.filename.decode('UTF-8')
    bam_file = pysam.AlignmentFile(bam_path, 'rb')

    # Beginning assembly
    if isinstance(merged_assembly, str):
        click.echo('Using the merged assembly as the reference insertion sequence...')
        right_assembly = merged_assembly
        left_assembly = merged_assembly
    else:
        click.echo('Using left and right flanks as separate reference sequecnes...')

    # Now get direct insertion sequence reads
    direct_inseq_reads = hndl_get_direct_inseq_reads(bam_file, contig, site, right_assembly, left_assembly,
                                                     output_prefix, outdir, site_name, out_dict)

    # Now get unique site ancestral reads
    direct_ancestral_reads = hndl_get_ancestral_inseq_reads(bam_file, contig, site, right_assembly, left_assembly,
                                                            output_prefix, outdir, site_name, out_dict)

    # Estimate the direct inseq frequency
    hndl_estimate_direct_inseq_freq(len(direct_inseq_reads), len(direct_ancestral_reads), out_dict)

    # Estimate the adjusted direct inseq frequency
    hndl_estimate_adj_direct_inseq_freq(max_softclip_length, average_read_length, len(direct_inseq_reads),
                                        len(direct_ancestral_reads), out_dict)

    del direct_ancestral_reads
    del direct_inseq_reads

    # Now get reads that flank the insertion site and map to the insertion
    flanking_inseq_reads = hndl_get_flanking_inseq_reads(bam_file, contig, site, flank_length,
                                                         right_assembly, left_assembly,
                                                         max_mapping_distance, output_prefix, outdir,
                                                         site_name, out_dict)

    # Now get the ancestral reads that flank the insertion...
    flanking_ancestral_reads = hndl_get_flanking_ancestral_reads(bam_file, contig, site, flank_length, right_assembly,
                                                                 left_assembly, max_mapping_distance,
                                                                 output_prefix, outdir, site_name, out_dict)

    # Quickly calculate the flanking inseq frequency
    hndl_estimate_flanking_inseq_freq(len(flanking_inseq_reads), len(flanking_ancestral_reads), out_dict)

    del flanking_inseq_reads
    del flanking_ancestral_reads

    return pd.DataFrame(out_dict)


def hndl_get_direct_inseq_reads(bam_file, contig, site, right_assembly, left_assembly, output_prefix, outdir, site_name, out_dict):
    click.echo('\tGetting IS overlapping reads directly at the insertion sequence site...')
    direct_inseq_reads = get_direct_inseq_reads(bam_file, contig, site, right_assembly, left_assembly)
    out_dict['direct_inseq_read_count'] = [len(direct_inseq_reads)]
    click.echo('\t\t%d IS overlapping reads directly at insertion site' % len(direct_inseq_reads))
    click.echo('\t\tWriting direct IS reads to file...')
    output.write_site_reads(direct_inseq_reads, bam_file, join(outdir, output_prefix + '.direct_inseq_bam'),
                            output_prefix, site_name + '.direct_inseq')
    click.echo()

    return direct_inseq_reads

def hndl_get_ancestral_inseq_reads(bam_file, contig, site, right_inseq_assembly, left_inseq_assembly, output_prefix,
                                   outdir, site_name, out_dict):
    click.echo('\tGetting ancestral reads directly overlapping site...')
    direct_ancestral_reads = get_direct_ancestral_reads(bam_file, contig, site, right_inseq_assembly,
                                                        left_inseq_assembly)
    out_dict['direct_ancestral_read_count'] = [len(direct_ancestral_reads)]
    click.echo('\t\t%d ancestral reads directly at insertion site' % len(direct_ancestral_reads))
    click.echo('\t\tWriting direct ancestral reads to file...')
    output.write_site_reads(direct_ancestral_reads, bam_file, join(outdir, output_prefix + '.direct_ancestral_bam'),
                            output_prefix, site_name + '.direct_ancestral')
    click.echo()

    return direct_ancestral_reads

def hndl_estimate_direct_inseq_freq(num_direct_inseq_reads, num_direct_ancestral_reads, out_dict):

    click.echo('\tCalculating the direct IS frequency...')
    click.echo('\t\tFirst caculating the raw direct IS frequency...')
    total_direct_reads = num_direct_inseq_reads + num_direct_ancestral_reads
    if total_direct_reads > 0:
        raw_direct_ineq_freq = num_direct_inseq_reads / total_direct_reads
    else:
        raw_direct_ineq_freq = np.nan
    out_dict['direct_inseq_freq'] = [raw_direct_ineq_freq]
    click.echo('\t\tThe raw direct IS frequency is %f\n' % raw_direct_ineq_freq)

def hndl_estimate_adj_direct_inseq_freq(max_softclip_length, average_read_length, num_direct_inseq_reads,
                                        num_direct_ancestral_reads, out_dict):
    click.echo('\t\tNow caculating the adjusted direct IS frequency...')
    click.echo('\t\tUsing a maximum softclip length of %d...' % max_softclip_length)
    direct_repeat_length = out_dict['direct_repeat_length'][0]
    click.echo('\t\tUsing direct repeat length of %d...' % direct_repeat_length)
    click.echo('\t\tUsing average read length of %d...' % average_read_length)
    click.echo('\t\tCalculating insertion sequence read count scaling factor...')
    inseq_scaling_factor = calculate_direct_inseq_scaling_factor(num_direct_inseq_reads, average_read_length,
                                                                 max_softclip_length)
    click.echo('\t\tThe scaling factor is %f...' % inseq_scaling_factor)
    click.echo('\t\tCalculating ancestral read count scaling factor...')
    ancestral_scaling_factor = calculate_direct_ancestral_scaling_factor(num_direct_ancestral_reads,
                                                                         average_read_length, direct_repeat_length)
    click.echo('\t\tThe scaling factor is %f...' % ancestral_scaling_factor)
    adjusted_inseq_read_count = num_direct_inseq_reads + inseq_scaling_factor * 4 * (
    average_read_length - max_softclip_length)
    adjusted_ancestral_read_count = 2 * num_direct_ancestral_reads + ancestral_scaling_factor * 4 * (
    direct_repeat_length + 1)
    adj_total_read_count = adjusted_inseq_read_count + adjusted_ancestral_read_count

    if adj_total_read_count > 0:
        adj_direct_inseq_freq = adjusted_inseq_read_count / adj_total_read_count
    else:
        adj_direct_inseq_freq = np.nan
    out_dict['adj_direct_inseq_freq'] = [adj_direct_inseq_freq]
    click.echo('\t\tThe adjusted direct IS frequency is %f' % adj_direct_inseq_freq)
    click.echo()

def hndl_get_flanking_inseq_reads(bam_file, contig, site, flank_length, right_assembly, left_assembly, max_mapping_distance,
                                  output_prefix, outdir, site_name, out_dict):
    click.echo('\tGetting reads that flank the IS site, and map to the relevant sequence...')
    flanking_inseq_reads = get_flanking_inseq_reads(bam_file, contig, site, flank_length, right_assembly,
                                                    left_assembly, max_mapping_distance)
    out_dict['flanking_inseq_read_count'] = [len(flanking_inseq_reads)]
    click.echo('\t\t%d reads flanking the insertion site' % len(flanking_inseq_reads))
    click.echo('\t\tWriting flanking IS reads to file...')
    output.write_site_reads(flanking_inseq_reads, bam_file, join(outdir, output_prefix + '.flanking_inseq_bam'),
                            output_prefix, site_name + '.flanking_inseq')
    click.echo()

    return flanking_inseq_reads

def hndl_get_flanking_ancestral_reads(bam_file, contig, site, flank_length, right_assembly, left_assembly, max_mapping_distance,
                                  output_prefix, outdir, site_name, out_dict):
    click.echo('\tGetting ancestral reads that flank the IS site...')
    flanking_ancestral_reads = get_flanking_ancestral_reads(bam_file, contig, site, flank_length,
                                                            right_assembly, left_assembly, max_mapping_distance)
    out_dict['flanking_ancestral_read_count'] = [len(flanking_ancestral_reads)]
    click.echo('\t\t%d ancestral reads flanking the insertion site' % len(flanking_ancestral_reads))
    click.echo('\t\tWriting flanking ancestral reads to file...')
    output.write_site_reads(flanking_ancestral_reads, bam_file,
                            join(outdir, output_prefix + '.flanking_ancestral_bam'), output_prefix,
                            site_name + '.flanking_ancestral')
    click.echo()

    return flanking_ancestral_reads

def hndl_estimate_flanking_inseq_freq(num_flanking_inseq_reads, num_flanking_ancestral_reads, out_dict):
    click.echo('\tCalculating the flanking IS frequency...')
    total_direct_reads = float(
        num_flanking_inseq_reads + num_flanking_ancestral_reads)
    if total_direct_reads > 0:
        flanking_inseq_freq = num_flanking_inseq_reads / total_direct_reads
    else:
        flanking_inseq_freq = np.nan
    out_dict['flanking_inseq_freq'] = [flanking_inseq_freq]
    click.echo('\t\tThe flanking IS frequency is %f' % flanking_inseq_freq)
    click.echo()