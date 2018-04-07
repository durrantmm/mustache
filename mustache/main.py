import sys
import click
import pysam
from mustache import bwa_tools, identify_candidate_sites, inseq_site_processor, misc, output
from os.path import join, basename, dirname, isdir
import os, sys
import pandas as pd
from Bio import SeqIO

@click.group()
def cli():
    pass

@click.group()
def align():
    pass

@click.group()
def find():
    pass

@click.group()
def callsites():
    pass

@click.command(name='paired')
@click.argument('fastq1', type=click.Path(exists=True))
@click.argument('fastq2', type=click.Path(exists=True))
@click.argument('genome', type=click.Path(exists=True))
@click.argument('out_bam')
@click.option('--threads', default=1, help="Specify the number of threads to use for the alignment.")
@click.option('--keep_tmp_files', is_flag=True, help="Add this flag if you want to keep intermediate temporary files")
def align_paired(fastq1, fastq2, genome, out_bam, threads, keep_tmp_files):
    click.echo("Performing alignment with bwa...")

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

    click.echo("Performing the alignment with BWA...")
    tmp_sam = join(dirname(out_bam), '.'.join(basename(out_bam).split('.')[:-1])+'.sam.tmp')
    aligned = bwa_tools.align_to_genome(fastq1, fastq2, genome, tmp_sam, threads)
    if not aligned:
        click.echo("Fatal error: Alignment seemed to have failed.")
        sys.exit()
    click.echo("Initial alignment completed successfully...\n")

    click.echo("Removing secondary alignments...")
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
    click.echo("SAM file successfully reformatted...\n")

    click.echo("Sorting the BAM file by chromosomal location...")
    sorted_bam = bwa_tools.samtools_sort_coordinate(tmp_formatted_bam, out_bam, delete_in_bam=not keep_tmp_files)
    if not sorted_bam:
        click.echo("Fatal error: Failed to sort the BAM file.")
        sys.exit()
    click.echo("BAM file successfully sorted...\n")

    click.echo("Index the sorted BAM file...")
    indexed = bwa_tools.samtools_index(out_bam)
    if not indexed:
        click.echo("Fatal error: Failed to index sorted BAM file")
        sys.exit()
    click.echo("BAM file successfully indexed...\n")


@click.command(name='single')
@click.argument('fastq', type=click.Path(exists=True))
@click.argument('genome', type=click.Path(exists=True))
@click.argument('out_bam')
@click.option('--threads', default=1, help="Specify the number of threads to use for the alignment.")
@click.option('--keep_tmp_files', is_flag=True, help="Add this flag if you want to keep intermediate temporary files")
def align_single(fastq, genome, out_bam, threads, keep_tmp_files):

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

    click.echo("Performing the alignment with bwa...")
    tmp_sam = join(dirname(out_bam), '.'.join(basename(out_bam).split('.')[:-1])+'.sam.tmp')
    aligned = bwa_tools.align_to_genome_single(fastq, genome, tmp_sam, threads)
    if not aligned:
        click.echo("Fatal error: Alignment seemed to have failed.")
        sys.exit()
    click.echo("Initial alignment completed successfully...\n")

    click.echo("Removing secondary alignments...")
    tmp_cleaned_bam = join(dirname(out_bam), '.'.join(basename(out_bam).split('.')[:-1]) + '.cleaned.bam.tmp')
    sorted_query_name = bwa_tools.samtools_remove_secondary_alignments(tmp_sam, tmp_cleaned_bam, delete_in_bam=not keep_tmp_files)
    if not sorted_query_name:
        click.echo("Fatal error: Failed to remove secondary alignments...")
        sys.exit()
    click.echo("Successfully removed secondary alignments...\n")

    click.echo("Sorting the BAM file by chromosomal location...")
    sorted_bam = bwa_tools.samtools_sort_coordinate(tmp_cleaned_bam, out_bam, delete_in_bam=not keep_tmp_files)
    if not sorted_bam:
        click.echo("Fatal error: Failed to sort the BAM file.")
        sys.exit()
    click.echo("BAM file successfully sorted...\n")

    click.echo("Index the sorted BAM file...")
    indexed = bwa_tools.samtools_index(out_bam)
    if not indexed:
        click.echo("Fatal error: Failed to index sorted BAM file")
        sys.exit()
    click.echo("BAM file successfully indexed...\n")


@click.command(name='paired')
@click.argument('bam_file', type=click.Path(exists=True))
@click.argument('genome_fasta', type=click.Path(exists=True))
@click.argument('output_prefix')
@click.option('--outdir', help="Indicate the directory where you want to write the output files. Default is the output given output prefix")
@click.option('--contig', help="Specify the contig you'd like to search.")
@click.option('--start', type=int, help="Specify the start point for your search.")
@click.option('--stop', type=int, help="Specify the stop point for your search")
@click.option('--min_softclip_length', default=4, help="Choose the minimum length for a site to be considered softclipped.")
@click.option('--min_softclip_count', default=4 , help="The minimum number of softclips to be considered a true signal.")
@click.option('--min_softclip_pair_distance', default=0, help="Choose the minimum distance between softclip sites in a pair.")
@click.option('--max_softclip_pair_distance', default=20, help="Choose the maximum distance between softclip sites in a pair.")
@click.option('--max_softclip_count_ratio_deviation', default=0.25, help="The maximum deviation from 0.5 allowable for softclip pairs.")
@click.option('--max_mapping_distance', default=3000, help="The maximum distance between two pairs to be considered concordantly mapped.")
@click.option('--flank_length', type=int, help="Specify flank length if desired. Otherwise, flank_length_percentile will be used to calculate flank length from the bam file.")
@click.option('--flank_length_percentile', default=0.999, help="Assuming normally distributed fragment size, choose the percentile used as length cutoff.")
@click.option('--assembly_kmer_size', default=21, help="Kmer size to use for velvet assemblies.")
@click.option('--min_inseq_read_count', default=4, help="After assembling a flank, minimum number of realigned reads to keep flank.")
@click.option('--threads', default=1, help="Specify the number of threads to use.")
def find_paired(bam_file, genome_fasta, output_prefix, outdir, contig, start, stop, min_softclip_length, min_softclip_count,
         min_softclip_pair_distance, max_softclip_pair_distance, max_softclip_count_ratio_deviation,
         max_mapping_distance, flank_length, flank_length_percentile, assembly_kmer_size, min_inseq_read_count, threads):

    if output_prefix.count('/') > 0:
        click.echo("Output prefix cannot be path to a different directory. Try using --outdir.")
        sys.exit()

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
    candidate_sites = identify_candidate_sites.get_softclipped_sites(bam_file, contig, start, stop, min_softclip_length,
                                                                     min_softclip_count,min_softclip_pair_distance,
                                                                     max_softclip_pair_distance)

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

    click.echo('Calculating the average read length from the BAM file...')
    avg_read_length = misc.calculate_average_read_length(bam_file)
    click.echo('Final average read length calculated as %d.\n' % avg_read_length )

    click.echo("Performing alignment-anchored assemblies of flanking sequencies...")
    assembled_flanks_df = None
    genome_name = None
    genome_seq = None

    for i in range(len(candidate_sites)):
        click.echo('\nProcessing site %d of %d' % (i+1, len(candidate_sites)))

        contig = candidate_sites.loc[i,'contig']
        site = candidate_sites.loc[i,'site']
        L_pass = candidate_sites.loc[i,'L_pass']
        R_pass = candidate_sites.loc[i, 'R_pass']

        # Retrieve the relevant contig from the fasta
        if genome_name != contig:
            click.echo("\tLoading relevant genome from fasta...")
            read_genome_fasta = SeqIO.parse(genome_fasta, format='fasta')

            for rec in read_genome_fasta:
                if rec.name == contig:
                    genome_name = rec.name
                    genome_seq = str(rec.seq)
                    break

        if L_pass:

            df = inseq_site_processor.process_candidate_site(bam_file, genome_seq, contig, site, 'L', flank_length,
                                                             max_mapping_distance, avg_read_length, assembly_kmer_size,
                                                             outdir, output_prefix)

            if not df is None:

                if assembled_flanks_df is None:
                    assembled_flanks_df = df
                else:
                    assembled_flanks_df = pd.concat([assembled_flanks_df, df])

        if R_pass:

            df = inseq_site_processor.process_candidate_site(bam_file, genome_seq, contig, site, 'R', flank_length,
                                                             max_mapping_distance, avg_read_length, assembly_kmer_size,
                                                             outdir, output_prefix)

            if not df is None:

                if assembled_flanks_df is None:
                    assembled_flanks_df = df
                else:
                    assembled_flanks_df = pd.concat([assembled_flanks_df, df])

    assembled_flanks_df.query("inseq_read_count >= @min_inseq_read_count", inplace=True)
    assembled_flanks_df = misc.keep_sites_with_nearby_mates2(assembled_flanks_df, min_softclip_pair_distance, max_softclip_pair_distance)
    assembled_flanks_df.reset_index(drop=True, inplace=True)

    click.echo()
    click.echo("Initial Assembled Flanks:")
    click.echo(assembled_flanks_df)

    site_pairs = inseq_site_processor.determine_site_pairs(assembled_flanks_df, bam_file, genome_fasta, min_softclip_pair_distance, max_softclip_pair_distance, outdir)

    priority_site_pairs = inseq_site_processor.prioritize_site_pairs(site_pairs)

    final_table = inseq_site_processor.create_final_output_table(assembled_flanks_df, priority_site_pairs)

    click.echo()
    click.echo("FINAL INSERTION SEQUENCE TABLE:")
    click.echo(final_table)

    output.write_final_dataframe(final_table, outdir, output_prefix)
    output.write_final_dataframe_to_fasta(final_table, outdir, output_prefix)



# Now set up the click arguments
align.add_command(align_paired)
align.add_command(align_single)

find.add_command(find_paired)

cli.add_command(align)
cli.add_command(find)


if __name__ == '__main__':
    cli()