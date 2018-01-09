from snakemake import shell
import pysam
from glob import glob
from os.path import isfile
import sys
import click


def genome_is_indexed(genome_path):

    amb = genome_path + '.amb'
    ann = genome_path + '.ann'
    bwt = genome_path + '.bwt'
    pac = genome_path + '.pac'
    sa = genome_path + '.sa'

    index_files = [amb, ann, bwt, pac, sa]

    indexed = True
    files = glob(genome_path+'*')
    for f in index_files:
        if f not in files:
            indexed = False

    return indexed


def index_genome(genome_path):
   shell('bwa index {genome_path}'.format(genome_path=genome_path))
   return genome_is_indexed(genome_path)


def align_to_genome(fastq1, fastq2, genome_path, out_sam, threads=1):
    command = "bwa mem -t {threads} {genome_path} {fastq1} {fastq2} > {out_sam}".format(
        genome_path=genome_path, fastq1=fastq1, fastq2=fastq2, out_sam=out_sam, threads=threads)
    click.echo("Executing command: %s" % command)
    shell(command)

    if isfile(out_sam):
        return True
    else:
        return False

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

def samtools_remove_secondary_alignments(in_bam, out_bam, delete_in_bam=False):
    shell('samtools view -b -h -F 0x900 {in_bam} > {out_bam}'.format(in_bam=in_bam, out_bam=out_bam))

    if delete_in_bam:
        shell('rm {in_bam}'.format(in_bam=in_bam))

    if isfile(out_bam):
        return True
    else:
        return False

def samtools_sort_name(in_bam, out_bam, delete_in_bam=False):
    shell('samtools sort -n {in_bam} > {out_bam}'.format(in_bam=in_bam, out_bam=out_bam))

    if delete_in_bam:
        shell('rm {in_bam}'.format(in_bam=in_bam))

    if isfile(out_bam):
        return True
    else:
        return False


def samtools_sort_coordinate(in_bam, out_bam, delete_in_bam=False):
    shell('samtools sort {in_bam} > {out_bam}'.format(in_bam=in_bam, out_bam=out_bam))

    if delete_in_bam:
        shell('rm {in_bam}'.format(in_bam=in_bam))

    if isfile(out_bam):
        return True
    else:
        return False


def samtools_index(in_bam):
    shell('samtools index {in_bam}'.format(in_bam=in_bam))

    if isfile(in_bam+'.bai'):
        return True
    else:
        return False