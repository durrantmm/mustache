import pysam
from os.path import join
from snakemake import shell
from bwa_tools import samtools_sort_coordinate, samtools_index
import sys
from os.path import isfile

def write_site_reads(reads, bam_file, outdir, output_prefix, suffix=None, sort_and_index=True):

    outfile_name = '.'.join([join(outdir, output_prefix), suffix, 'bam'])

    if sort_and_index:
        tmp_outfile_name = outfile_name+'.tmp'
        tmp_out_bam = pysam.AlignmentFile(tmp_outfile_name, 'wb', template=bam_file)
        for r in reads:
            tmp_out_bam.write(r)
        tmp_out_bam.close()

        samtools_sort_coordinate(tmp_outfile_name, outfile_name, delete_in_bam=True)
        samtools_index(outfile_name)
    else:
        out_bam = pysam.AlignmentFile(outfile_name, 'wb', template=bam_file)
        for r in reads:
            out_bam.write(r)
        out_bam.close()

    if isfile(outfile_name):
        return True
    else:
        return False

def write_final_dataframe(df, outdir, output_prefix):
    outfile_name = join(outdir, output_prefix) + '.stats.csv'
    df.to_csv(outfile_name, index=False)