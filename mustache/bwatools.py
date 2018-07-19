import sys
from glob import glob
from snakemake import shell
import pygogo as gogo
from os.path import isfile

verbose=False
logger = gogo.Gogo(__name__, verbose=verbose).logger

def index_genome(genome_path, silence=True):
    if silence:
        shell('bwa index {genome_path} 2> /dev/null;'.format(genome_path=genome_path))
    else:
        shell('bwa index {genome_path}'.format(genome_path=genome_path))

    return genome_is_indexed(genome_path)

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

def align_to_genome_pe(fastq1, fastq2, genome_path, out_sam, threads=1, verbose=False):
    if verbose:
        command = "bwa mem -t {threads} {genome_path} {fastq1} {fastq2} > {out_sam}; ".format(
            genome_path=genome_path, fastq1=fastq1, fastq2=fastq2, out_sam=out_sam, threads=threads)
    else:
        command = "bwa mem -t {threads} {genome_path} {fastq1} {fastq2} 2> /dev/null 1> {out_sam}; ".format(
            genome_path=genome_path, fastq1=fastq1, fastq2=fastq2, out_sam=out_sam, threads=threads)

    logger.debug("Executing command: %s" % command)
    shell(command)

    if isfile(out_sam):
        return True
    else:
        return False

def align_to_genome_se(fastq1, genome_path, out_sam, threads=1, verbose=False):
    if verbose:
        command = "bwa mem -t {threads} {genome_path} {fastq1} > {out_sam}; ".format(
            genome_path=genome_path, fastq1=fastq1, out_sam=out_sam, threads=threads)
    else:
        command = "bwa mem -t {threads} {genome_path} {fastq1} 2> /dev/null 1> {out_sam}; ".format(
            genome_path=genome_path, fastq1=fastq1, out_sam=out_sam, threads=threads)

    logger.debug("Executing command: %s" % command)
    shell(command)

    if isfile(out_sam):
        return True
    else:
        return False