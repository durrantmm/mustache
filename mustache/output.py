import pysam
from mustache.bwa_tools import samtools_sort_coordinate, samtools_index
import sys
from snakemake import shell
from os.path import isfile
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


def write_final_dataframe(df, outdir, output_prefix, suffix='.insertion_seqs.tsv'):
    outfile_name = join(outdir, output_prefix) + suffix
    df.to_csv(outfile_name, index=False, sep='\t')

def write_final_dataframe_to_fasta(df, outdir, output_prefix, suffix='.insertion_seqs.fasta'):

    outfile_name = join(outdir, output_prefix) + suffix

    out_sequences = []
    for index, row in df.iterrows():

        contig = row['contig']
        left_site = row['left_site']
        right_site = row['right_site']
        left_read_counts = row['left_site_read_count']
        right_read_counts = row['right_site_read_count']
        assembly_length = row['assembly_length']
        assembly = row['assembly']
        partner_site = row['partner_site']

        if row['orientation'] == 'M':
            name = '_'.join(map(str, [contig, int(left_site), int(right_site), 'MERGED', 'readcount',
                                      int(left_read_counts), int(right_read_counts), 'length', int(assembly_length)]))
        else:
            if row['orientation'] == 'L':
                name = '_'.join(map(str, [contig, int(left_site), 'LEFT_FLANK', 'partner', int(partner_site),
                                          'readcount', int(left_read_counts), 'length', int(assembly_length)]))
            else:
                name = '_'.join(map(str, [contig, int(right_site), 'RIGHT_FLANK', 'partner', int(partner_site),
                                          'readcount', int(right_read_counts), 'length', int(assembly_length)]))

        record = SeqRecord(Seq(assembly, IUPAC.IUPACAmbiguousDNA),
                           id=name, description=name)

        out_sequences.append(record)

    with open(outfile_name, 'w') as output_handle:
        SeqIO.write(out_sequences, output_handle, 'fasta')


def write_final_merged_dataframe_to_fasta(df, outdir, output_prefix, suffix='.merged_insertion_seqs.fasta'):

    outfile_name = join(outdir, output_prefix) + suffix

    out_sequences = []
    for index, row in df.iterrows():

        contig = row['contig']
        left_site = row['left_site']
        right_site = row['right_site']
        assembly_length = row['assembly_length']
        assembly = row['assembly']
        partner_site = row['partner_site']

        if row['orientation'] == 'M':
            name = '_'.join(map(str, [contig, int(left_site), int(right_site), 'MERGED', 'length', int(assembly_length)]))
        else:
            if row['orientation'] == 'L':
                name = '_'.join(map(str, [contig, int(left_site), 'LEFT_FLANK', 'partner', int(partner_site),
                                          'length', int(assembly_length)]))
            else:
                name = '_'.join(map(str, [contig, int(right_site), 'RIGHT_FLANK', 'partner', int(partner_site),
                                          'length', int(assembly_length)]))

        record = SeqRecord(Seq(assembly, IUPAC.IUPACAmbiguousDNA),
                           id=name, description=name)

        out_sequences.append(record)

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
    command_merge = 'samtools merge -r -f {out_bam} {in_bams}'.format(out_bam=out_bam, in_bams=in_bams)
    shell(command_merge)
    command_index = 'samtools index {out_bam}'.format(out_bam=out_bam)
    shell(command_index)

def dataframe_to_stdout(df):
    df.to_csv(sys.stdout, index=False, sep='\t')


def write_reads_to_fasta(reads, fasta, name_index=0, seq_index=1):

    with open(fasta, 'w') as outfasta:
        for r in reads:
            outfasta.write(">" + r[name_index] + '\n')
            outfasta.write(r[seq_index] + '\n')