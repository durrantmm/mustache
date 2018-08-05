import warnings
warnings.filterwarnings("ignore")
import sys
from snakemake import shell
from Bio import SeqIO
import pandas as pd
from mustache.config import FRAGGENESCAN
#conda install -c bioconda fraggenescan 

def run_fraggenescan(fasta, outfile, threads=2):
    fraggenescan_command = '{fraggenescan} -genome={fasta} -out={outfile} -complete=0 -train=complete -thread={threads} 1> /dev/null'.format(
        fraggenescan=FRAGGENESCAN, fasta=fasta, outfile=outfile, threads=threads
    )
    shell(fraggenescan_command)

def count_fraggenescan_seqs(fasta):
    count = 0
    for rec in SeqIO.parse(fasta, 'fasta'):
        count += 1
    return count