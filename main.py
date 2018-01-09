import sys
import click
import pysam
import bwa_tools, identify_candidate_sites, inseq_site_processor, misc
from os.path import join, basename, dirname, isdir, isfile
import os, sys
import pandas as pd
import output
from multiprocessing import Pool

@click.group()
def cli():
    pass

@click.command()
@click.argument('fastq1', type=click.Path(exists=True))
@click.argument('fastq2', type=click.Path(exists=True))
@click.argument('genome', type=click.Path(exists=True))
@click.argument('out_bam')
@click.option('--threads', default=1, help="Specify the number of threads to use for the alignment.")
@click.option('--keep_tmp_files', is_flag=True, help="Add this flag if you want to keep intermediate temporary files")
def align(fastq1, fastq2, genome, out_bam, threads, keep_tmp_files):
    click.echo("Performing alignment with bowtie2...")

    click.echo("Checking if genome is indexed...")
    if bwa_tools.genome_is_indexed(genome):
        click.echo("Genome is already indexed, skipping...")
    else:
        click.echo("Indexing Genome...")
        genome_indexed = bwa_tools.index_genome(genome)
        if not genome_indexed:
            click.echo("Fatal error: Genome could not be indexed. Is it in FASTA format?")
            sys.exit()
        click.echo("The genome is properly indexed...\n")

    click.echo("Performing the alignment with bowtie2...")
    tmp_sam = join(dirname(out_bam), '.'.join(basename(out_bam).split('.')[:-1])+'.sam.tmp')
    aligned = bwa_tools.align_to_genome(fastq1, fastq2, genome, tmp_sam, threads)
    if not aligned:
        click.echo("Fatal error: Alignment seemed to have failed.")
        sys.exit()
    click.echo("Initial alignment completed successfully...\n")

    click.echo("Removing secondary alignmnents...")
    tmp_cleaned_bam = join(dirname(out_bam), '.'.join(basename(out_bam).split('.')[:-1]) + '.cleaned.bam.tmp')
    sorted_query_name = bwa_tools.samtools_remove_secondary_alignments(tmp_sam, tmp_cleaned_bam, delete_in_bam=not keep_tmp_files)
    if not sorted_query_name:
        click.echo("Fatal error: Failed to remove secondary alignments...")
        sys.exit()
    click.echo("Successfully removed secondary alignments...\n")

    click.echo("Formatting bam file for use by mustache...")
    tmp_formatted_bam = join(dirname(out_bam), '.'.join(basename(out_bam).split('.')[:-1]) + '.formatted.bam.tmp')
    reformatted = bwa_tools.format_for_mustache(tmp_cleaned_bam, tmp_formatted_bam, delete_in_sam=not keep_tmp_files)
    if not reformatted:
        click.echo("Fatal error: SAM file reformatting failed.")
        sys.exit()
    click.echo("SAM file successfully reformated...\n")

    click.echo("Sorting the BAM file by chromosomal location...")
    sorted = bwa_tools.samtools_sort_coordinate(tmp_formatted_bam, out_bam, delete_in_bam=not keep_tmp_files)
    if not sorted:
        click.echo("Fatal error: Failed to sort the BAM file.")
        sys.exit()
    click.echo("BAM file successfully sorted...\n")

    click.echo("Index the sorted BAM file...")
    indexed = bwa_tools.samtools_index(out_bam)
    if not indexed:
        click.echo("Fatal error: Failed to index sorted BAM file")
        sys.exit()
    click.echo("BAM file successfully indexed...\n")


@click.command()
@click.argument('bam_file', type=click.Path(exists=True))
@click.argument('output_prefix')
@click.option('--outdir', help="Indicate the directory where you want to write the output files. Default is the output given output prefix")
@click.option('--contig', help="Specify the contig you'd like to search.")
@click.option('--start', type=int, help="Specify the start point for your search.")
@click.option('--stop', type=int, help="Specify the stop point for your search")
@click.option('--min_softclip_length', default=5, help="Choose the minimum length for a site to be considered softclipped.")
@click.option('--min_softclip_count', default=5, help="The minimum number of softclips to be considered a true signal.")
@click.option('--min_softclip_pair_distance', default=2, help="Choose the minimum distance between softclip sites in a pair.")
@click.option('--max_softclip_pair_distance', default=15, help="Choose the maximum distance between softclip sites in a pair.")
@click.option('--max_softclip_count_ratio_deviation', default=0.25, help="The maximum deviation from 0.5 allowable for softclip pairs.")
@click.option('--max_mapping_distance', default=3000, help="The maximum distance between two pairs to be considered concordantly mapped.")
@click.option('--flank_length', type=int, help="Specify flank length if desired. Otherwise, flank_length_percentile will be used to calculate flank length from the bam file.")
@click.option('--flank_length_percentile', default=0.999, help="Assuming normally distributed fragment size, choose the percentile used as length cutoff.")
@click.option('--threads', default=1, help="Specify the number of threads to use.")
def find(bam_file, output_prefix, outdir, contig, start, stop, min_softclip_length, min_softclip_count,
         min_softclip_pair_distance, max_softclip_pair_distance, max_softclip_count_ratio_deviation,
         max_mapping_distance, flank_length, flank_length_percentile, threads):

    if not outdir:
        outdir = output_prefix
    else:
        outdir = outdir

    click.echo("Making output directory if it doesn't exist...")
    os.makedirs(outdir, exist_ok=True)
    if not isdir(outdir):
        click.echo("Fatal error: Could not create specified directory.")
        sys.exit()
    click.echo("Successfully created the specified output directory...\n")

    click.echo("Opening specified BAM file...")
    bam_file = pysam.AlignmentFile(bam_file, 'rb')
    click.echo("BAM file successfully opened...\n")

    click.echo("Identifying the candidate insertion sites...")
    candidate_sites = identify_candidate_sites.run(bam_file, contig, start, stop, min_softclip_length, min_softclip_count,
                                                   min_softclip_pair_distance, max_softclip_pair_distance,
                                                   max_softclip_count_ratio_deviation)
    if len(candidate_sites) == 0:
        click.echo("No candidate sites identified. Try different parameters.")
        sys.exit()
    click.echo("A total of %d candidate sites have been identified...\n" % len(candidate_sites))
    click.echo(candidate_sites)
    click.echo()

    if not flank_length:
        click.echo('Calculating the maximum inseq flank length based on observed fragment length...')
        flank_length = misc.calculate_flank_length(bam_file, max_mapping_distance, flank_length_percentile)
        click.echo('Final flank length calculated as %d.' % flank_length)
        click.echo('Assuming fragment length to be normally distributed, this should include {perc}% of all read pairs...\n'.format(perc = 100*flank_length_percentile))
    else:
        click.echo('Using specified flank length of %d...\n' % flank_length)


    final_df = None
    for i in range(len(candidate_sites)):
        click.echo('Processing site %d of %d' % (i+1, len(candidate_sites)))
        contig = candidate_sites.loc[i,'contig']
        site = (candidate_sites.loc[i,'left_site'], candidate_sites.loc[i,'right_site'])
        df = inseq_site_processor.process_candidate_site(bam_file, contig, site, flank_length,
                                                         max_mapping_distance, outdir, output_prefix)
        if final_df is None:
            final_df = df
        else:
            final_df = pd.concat([final_df, df])

    output.write_final_dataframe(final_df, outdir, output_prefix)





cli.add_command(align)
cli.add_command(find)


