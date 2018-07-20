import sys
from snakemake import shell
import pandas as pd

#conda install -c bioconda fraggenescan 

def run_fraggenescan(fasta, outfile, threads=2):
    fraggenescan_command = 'run_FragGeneScan.pl -genome={fasta} -out={outfile} -complete=0 -train=complete -thread={threads}'.format(
        fasta=fasta, outfile=outfile, threads=threads
    )
    shell(fraggenescan_command)
