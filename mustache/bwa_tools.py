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


def index_genome(genome_path, silence=False):
    if silence:
        shell('bwa index {genome_path} 2> /dev/null;'.format(genome_path=genome_path))
    else:
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


def align_to_genome_single(fastq, genome_path, out_sam, additional_flags='', threads=1):
    command = "bwa mem -t {threads} {additional_flags} {genome_path} {fastq} > {out_sam}".format(
        genome_path=genome_path, fastq=fastq,out_sam=out_sam, threads=threads, additional_flags=additional_flags)
    click.echo("Executing command: %s" % command)
    shell(command)

    if isfile(out_sam):
        return True
    else:
        return False

def bwa_aln_to_genome_single(fastq, genome_path, out_bam, additional_flags='', threads=1, silence=False, qfilter=True):
    out_sai = out_bam+'.sai'

    if silence:
        if qfilter:
            command = "bwa aln -t {threads} {additional_flags} {genome_path} {fastq} 2> /dev/null 1> {out_sai}; " \
                      "bwa samse {genome_path} {out_sai} {fastq} 2> /dev/null 1> {out_bam}.tmp; " \
                      "rm {out_sai}; " \
                      "samtools sort {out_bam}.tmp 2> /dev/null 1> {out_bam}; " \
                      "samtools index {out_bam} 2> /dev/null; " \
                      "samtools view -h -q 1 {out_bam}.tmp 2> /dev/null 1> {out_bam}.tmp.tmp; " \
                      "samtools sort {out_bam}.tmp.tmp 2> /dev/null 1> {out_bam}; " \
                      "samtools index {out_bam}; " \
                      "rm -f {out_bam}.tmp {out_bam}.tmp.tmp".format(
                genome_path=genome_path, fastq=fastq, out_bam=out_bam, out_sai=out_sai, threads=threads, additional_flags=additional_flags)
        else:
            command = "bwa aln -t {threads} {additional_flags} {genome_path} {fastq} 2> /dev/null 1> {out_sai}; " \
                      "bwa samse {genome_path} {out_sai} {fastq} 2> /dev/null 1> {out_bam}.tmp; " \
                      "rm {out_sai}; " \
                      "samtools sort {out_bam}.tmp 2> /dev/null 1> {out_bam}; " \
                      "samtools index {out_bam} 2> /dev/null; " \
                      "samtools view -h {out_bam}.tmp 2> /dev/null 1> {out_bam}.tmp.tmp; " \
                      "samtools sort {out_bam}.tmp.tmp 2> /dev/null 1> {out_bam}; " \
                      "samtools index {out_bam}; " \
                      "rm -f {out_bam}.tmp {out_bam}.tmp.tmp".format(
                genome_path=genome_path, fastq=fastq, out_bam=out_bam, out_sai=out_sai, threads=threads,
                additional_flags=additional_flags)

    else:
        if qfilter:
            command = "bwa aln -t {threads} {additional_flags} {genome_path} {fastq} > {out_sai}; " \
                      "bwa samse {genome_path} {out_sai} {fastq} > {out_bam}.tmp; " \
                      "rm {out_sai}; " \
                      "samtools view -h -q 1 {out_bam}.tmp > {out_bam}.tmp.tmp; " \
                      "samtools sort {out_bam}.tmp.tmp > {out_bam}; " \
                      "samtools index {out_bam}; " \
                      "rm -f {out_bam}.tmp {out_bam}.tmp.tmp".format(
                genome_path=genome_path, fastq=fastq, out_bam=out_bam, out_sai=out_sai, threads=threads,
                additional_flags=additional_flags)
        else:
            command = "bwa aln -t {threads} {additional_flags} {genome_path} {fastq} > {out_sai}; " \
                      "bwa samse {genome_path} {out_sai} {fastq} > {out_bam}.tmp; " \
                      "rm {out_sai}; " \
                      "samtools view -h {out_bam}.tmp > {out_bam}.tmp.tmp; " \
                      "samtools sort {out_bam}.tmp.tmp > {out_bam}; " \
                      "samtools index {out_bam}; " \
                      "rm -f {out_bam}.tmp {out_bam}.tmp.tmp".format(
                genome_path=genome_path, fastq=fastq, out_bam=out_bam, out_sai=out_sai, threads=threads,
                additional_flags=additional_flags)

    if not silence:
        click.echo("Executing command: %s" % command)

    shell(command, read=silence)

    if isfile(out_bam):
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