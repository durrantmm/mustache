import sys
import click
from mustache import bwatools, samtools
from os.path import dirname, basename, join, isfile
import pysam
from snakemake import shell
import pygogo as gogo

verbose=True
logger = gogo.Gogo(__name__, verbose=verbose).logger


def read_sam_pairs(samfile):
    while True:
        yield(next(samfile), next(samfile))

def format_for_mustache(in_sam, out_bam, read_format='rb', delete_in_sam=False):

    samfile = pysam.AlignmentFile(in_sam, read_format)
    outbam= pysam.AlignmentFile(out_bam, "wb", template=samfile)
    sam_pairs = read_sam_pairs(samfile)

    for p1, p2 in sam_pairs:

        p1_qual = p1.tostring(samfile).split('\t')[10]
        p2_qual = p2.tostring(samfile).split('\t')[10]

        if p1.query_name != p2.query_name:

            click.echo("Fatal Error: Mismatched reads in bwa_tools.format()")
            sys.exit()

        if not p1.is_unmapped and p2.is_unmapped:
            p1.set_tag('MT', p2.query_sequence)
            p1.set_tag('MQ', p2_qual)
            outbam.write(p1)
        elif p1.is_unmapped and not p2.is_unmapped:
            p2.set_tag('MT', p1.query_sequence)
            p2.set_tag('MQ', p1_qual)
            outbam.write(p2)
        elif not p1.is_unmapped and not p2.is_unmapped:
            p1.set_tag('MT', p2.query_sequence)
            p1.set_tag('MQ', p2_qual)
            p2.set_tag('MT', p1.query_sequence)
            p2.set_tag('MQ', p1_qual)
            outbam.write(p1)
            outbam.write(p2)

    samfile.close()
    outbam.close()
    if delete_in_sam:
        shell('rm {in_sam}'.format(in_sam=in_sam))

    if isfile(out_bam):
        return out_bam
    else:
        return None


@click.command()
@click.argument('genome', type=click.Path(exists=True))
@click.argument('out_bam')
@click.argument('fastq1',  type=click.Path(exists=True))
@click.argument('fastq2', required=False, type=click.Path(exists=True))
@click.option('--threads', default=1, help="Specify the number of threads to use for the alignment.")
@click.option('--keep_tmp_files', is_flag=True, help="Add this flag if you want to keep intermediate temporary files")
def alignbwa(genome, out_bam, fastq1, fastq2, threads, keep_tmp_files):
    logger.info("Performing alignment with bwa...\n")

    logger.info("Checking if genome is indexed...")
    if bwatools.genome_is_indexed(genome):
        logger.info("Genome is already indexed, skipping...\n")
    else:
        logger.info("Indexing Genome...")
        genome_indexed = bwatools.index_genome(genome)
        if not genome_indexed:
            logger.error("Fatal error: Genome could not be indexed. Is it in FASTA format?")
            sys.exit()
        logger.info("The genome is properly indexed...\n")

    logger.info("Performing the alignment with BWA...")

    tmp_sam = join(dirname(out_bam), '.'.join(basename(out_bam).split('.')[:-1]) + '.bam.tmp')

    if type(fastq2) is type(None):
        aligned = bwatools.align_to_genome_se(fastq1, genome, tmp_sam, threads)
    else:
        aligned = bwatools.align_to_genome_pe(fastq1, fastq2, genome, tmp_sam, threads)


    if not aligned:
        logger.error("Fatal error: Alignment seemed to have failed.")
        sys.exit()

    logger.info("Initial alignment completed successfully...\n")

    logger.info("Removing secondary alignments...")
    tmp_cleaned_bam = join(dirname(out_bam), '.'.join(basename(out_bam).split('.')[:-1]) + '.cleaned.bam.tmp')
    sorted_query_name = samtools.remove_secondary_alignments(tmp_sam, tmp_cleaned_bam, delete_in_bam=not keep_tmp_files)
    if not sorted_query_name:
        logger.info("Fatal error: Failed to remove secondary alignments...")
        sys.exit()
    logger.info("Successfully removed secondary alignments...\n")

    if type(fastq2) is type(None):
        logger.info("Since file is single-ended, no need to format for mustache alignment-anchored assembly...\n")
        tmp_formatted_bam = tmp_cleaned_bam
    else:
        logger.info("Formatting bam file for use by mustache...")
        tmp_formatted_bam = join(dirname(out_bam), '.'.join(basename(out_bam).split('.')[:-1]) + '.formatted.bam.tmp')
        reformatted = format_for_mustache(tmp_cleaned_bam, tmp_formatted_bam, delete_in_sam=not keep_tmp_files)
        if not reformatted:
            logger.error("Fatal error: SAM file reformatting failed.")
            sys.exit()
        logger.info("SAM file successfully reformatted...\n")


    logger.info("Sorting the BAM file by chromosomal location...")
    sorted_bam = samtools.sort_coordinate(tmp_formatted_bam, out_bam, delete_in_bam=not keep_tmp_files)
    if not sorted_bam:
        logger.error("Fatal error: Failed to sort the BAM file.")
        sys.exit()
    logger.info("BAM file successfully sorted...\n")

    logger.info("Index the sorted BAM file...")
    indexed = samtools.index(out_bam)
    if not indexed:
        logger.info("Fatal error: Failed to index sorted BAM file")
        sys.exit()
    logger.info("BAM file successfully indexed...\n")


if __name__ == '__main__':
    alignbwa()