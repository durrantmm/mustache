import sys
from snakemake import shell
import os
from os.path import join
import subprocess
from Bio import SeqIO
import pysam
import time, timeit
from random import randint
from mustache.bwatools import index_genome, align_to_genome_se
from mustache.misc import revcomp

class MinimusAssembler:

    full_outprefix = None

    reads_path = None
    reads = None
    read_names = None
    afg_path = None
    out_fasta = None
    total_aligned_reads = None

    bam_path = None

    def __init__(self, reads, outdir='/tmp', outprefix='mustache.minimus.'+str(randint(0,1e100)), read_names=None):

        self.reads = reads
        self.read_names = read_names

        self.full_outprefix = join(outdir, outprefix)
        self.reads_path = join(outdir, outprefix + '.reads.fasta')
        self.align_seq_fasta_path = join(outdir, outprefix + '.align_seq.fasta')
        self.align_sam_path = join(outdir, outprefix + '.align_seq.sam')
        self.afg_path = join(outdir, outprefix + '.afg')
        self.bam_path = join(outdir, outprefix + '.bam')

        os.makedirs(outdir, exist_ok=True)

    def delete_files(self):
        shell('rm -rf {reads} {afg} {outprefix}.bnk {outprefix}.fasta'.format(outprefix=self.full_outprefix,
                                                                              afg=self.afg_path, reads=self.reads_path))

    def delete_afg_bank(self):
        shell('rm -rf {afg} {outprefix}.bnk'.format(outprefix=self.full_outprefix, afg=self.afg_path,
                                                    reads=self.reads_path))

    def assemble(self, min_overlap_length=10, min_overlap_count=2, stranded=False):
        self.write_reads_as_fasta()
        self.delete_afg_bank()
        shell("toAmos -s {reads_path} -o {afg_path} 2> /dev/null;".format(reads_path=self.reads_path,
                                                                          afg_path=self.afg_path), read=True)
        shell("bank-transact -f -z -b {out_prefix}.bnk -m {out_prefix}.afg 2> /dev/null;".format(
            out_prefix=self.full_outprefix), read=True)
        if stranded:
            shell(
                "hash-overlap -s -o {min_olen} -B {out_prefix}.bnk 2> /dev/null;".format(out_prefix=self.full_outprefix,
                                                                                         min_olen=min_overlap_length),
                read=True)
        else:
            shell("hash-overlap -o {min_olen} -B {out_prefix}.bnk 2> /dev/null;".format(out_prefix=self.full_outprefix,
                                                                                        min_olen=min_overlap_length),
                  read=True)
        shell("tigger -b {out_prefix}.bnk 2> /dev/null;".format(out_prefix=self.full_outprefix), read=True)
        shell(
            "make-consensus -o {min_ocount} -B -b {out_prefix}.bnk 2> /dev/null;".format(out_prefix=self.full_outprefix,
                                                                                         min_ocount=min_overlap_count),
            read=True)
        shell("bank2fasta -d -b {out_prefix}.bnk > {out_prefix}.fasta 2> /dev/null;".format(
            out_prefix=self.full_outprefix), read=True)

        self.out_fasta = join(self.full_outprefix + '.fasta')
        return self.out_fasta

    def write_reads_as_fasta(self):
        with open(self.reads_path, 'w') as out:
            for i in range(len(self.reads)):
                if self.read_names:
                    out.write('>%s' % self.read_names[i] + '\n' + self.reads[i] + '\n')
                else:
                    out.write('>Read%d' % i + '\n' + self.reads[i] + '\n')

    def get_total_aligned_reads(self):
        if self.total_aligned_reads:
            return self.total_aligned_reads

        total = 0
        for contig in SeqIO.parse(self.out_fasta, 'fasta'):
            nreads = int(contig.description.split()[-2].split('=')[-1])
            total += nreads
        self.total_aligned_reads = total
        return self.total_aligned_reads

    def get_assembled_sequences(self):
        for contig in SeqIO.parse(self.out_fasta, 'fasta'):
            yield contig


    def count_assembled_seqs(self):
        n_assembled_seqs = 0
        for contig in SeqIO.parse(self.out_fasta, 'fasta'):
            n_assembled_seqs += 1
        return n_assembled_seqs

    def write_seq_to_fasta(self, seq):
        with open(self.align_seq_fasta_path, 'w') as out:
           out.write('>align_seq\n' + seq)

    def align_seq_to_assembly(self, seq):
        self.write_seq_to_fasta(seq)
        index_genome(self.out_fasta)
        align_to_genome_se(self.align_seq_fasta_path, self.out_fasta, self.align_sam_path)

    def retrieve_extended_sequence(self, orient):
        outsam = pysam.AlignmentFile(self.align_sam_path, 'r')

        extension = None
        for read in outsam:
            if read.is_unmapped:
                return None

            query_length = len(read.query_sequence)

            current_contig = None
            for contig in self.get_assembled_sequences():
                if contig.id == read.reference_name:
                    current_contig = contig

            if orient == 'R' and (not read.is_reverse):
                if len(current_contig.seq) - read.reference_start <= query_length:
                    return None
                extension = current_contig.seq[read.reference_start:]
            elif orient == 'R' and read.is_reverse:
                if read.reference_start <= query_length:
                    return None
                extension = revcomp(current_contig.seq[:read.reference_end])
            elif orient == 'L' and (not read.is_reverse):
                if read.reference_start <= query_length:
                    return None
                extension = current_contig.seq[:read.reference_end]
            elif orient == 'L' and read.is_reverse:
                if len(current_contig.seq) - read.reference_start <= query_length:
                    return None
                extension = revcomp(current_contig.seq[read.reference_start:])

        return extension


