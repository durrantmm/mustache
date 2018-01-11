import pysam
from os.path import join
from bwa_tools import samtools_sort_coordinate, samtools_index
import sys
from snakemake import shell
from os.path import isfile
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import os
from os.path import join


def write_site_reads(reads, bam_file, outdir, output_prefix, suffix=None, sort_and_index=True):

    os.makedirs(outdir, exist_ok=True)

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


def stats_dataframe_to_fasta(df, outdir, output_prefix):

    outfile_name = join(outdir, output_prefix) + '.insertion_seqs.fna'
    out_sequences = []

    for row in range(len(df)):

        if np.isnan(df.iloc[row,]['merged_assembly_length']):

            contig = df.iloc[row,]['contig']
            left_site = df.iloc[row,]['left_site']
            right_site = df.iloc[row,]['right_site']

            right_assembly = df.iloc[row,]['right_assembly']
            left_assembly = df.iloc[row,]['left_assembly']

            right_assembly_name = '_'.join(map(str, [contig, left_site, right_site, 'right']))
            left_assembly_name = '_'.join(map(str, [contig, left_site, right_site, 'left']))

            right_assembly_record = SeqRecord(Seq(right_assembly, IUPAC.IUPACUnambiguousDNA),
                                              id=right_assembly_name, description=right_assembly_name)
            left_assembly_record = SeqRecord(Seq(left_assembly, IUPAC.IUPACUnambiguousDNA),
                                             id=left_assembly_name, description=left_assembly_name)

            out_sequences.append(right_assembly_record)
            out_sequences.append(left_assembly_record)

        elif df.iloc[row,]['merged_assembly_length'] != np.nan:

            contig = df.iloc[row,]['contig']
            left_site = df.iloc[row,]['left_site']
            right_site = df.iloc[row,]['right_site']

            merged_assembly = df.iloc[row,]['merged_assembly']

            merged_assembly_name = '_'.join(map(str, [contig, left_site, right_site, 'merged']))

            merged_assembly_record = SeqRecord(Seq(merged_assembly, IUPAC.IUPACUnambiguousDNA),
                                               id=merged_assembly_name, description=merged_assembly_name)

            out_sequences.append(merged_assembly_record)

    with open(outfile_name, 'w') as output_handle:
        SeqIO.write(out_sequences, output_handle, 'fasta')


def stats_dataframe_to_gff3(df, outdir, output_prefix):

    outfile_name = join(outdir, output_prefix) + '.sites.gff3'

    with open(outfile_name, 'w') as out:

        for row in range(len(df)):

            contig = df.iloc[row,]['contig']
            left_site = df.iloc[row,]['left_site']
            right_site = df.iloc[row,]['right_site']

            right_assembly = df.iloc[row,]['right_assembly']
            left_assembly = df.iloc[row,]['left_assembly']
            merged_assembly = df.iloc[row,]['merged_assembly']

            name = 'insertion_site-' + '_'.join(map(str, [contig, left_site, right_site]))

            out.write('\t'.join(map(str, [contig, 'mustache', 'insertion_site', left_site+2, right_site, name, '.',
                                          '.', '.', 'RIGHT_ASSEMBLY={right};LEFT_ASSEMBLY={left};MERGED_ASSEMBLY={merged}\n'.format(
                                        left=left_assembly, right=right_assembly, merged=merged_assembly)])))


def write_sequence(sequence, name, out_path):
    record = SeqRecord(Seq(sequence, IUPAC.IUPACUnambiguousDNA), id=name, description=name)

    with open(out_path, "w") as output_handle:
        SeqIO.write([record], output_handle, "fasta")


def merge_bams_in_directory(direc, out_bam):
    in_bams = join(direc, '*.bam')
    command_merge = 'samtools merge -f {out_bam} {in_bams}'.format(out_bam=out_bam, in_bams=in_bams)
    shell(command_merge)
    command_index = 'samtools index {out_bam}'.format(out_bam=out_bam)
    shell(command_index)
