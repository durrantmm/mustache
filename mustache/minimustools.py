import warnings
warnings.filterwarnings("ignore")
import sys
from snakemake import shell
import os
from os.path import join, isfile
from Bio import SeqIO
import pysam
from mustache.bwatools import index_genome, align_to_genome_se
from mustache.sctools import left_softclipped_sequence_strict, right_softclipped_sequence_strict
from mustache.misc import revcomp
from mustache.pysamtools import query_qualities_to_phred
from random import randint, sample


class MinimusAssembler:

    outdir = None
    full_outprefix = None

    reads_path = None
    reads = None
    read_names = None
    afg_path = None
    out_fasta = None
    total_aligned_reads = None

    bam_path = None

    def __init__(self, reads, quals=None, outdir=None, outprefix=None, read_names=None):

        self.outdir = outdir
        self.outprefix = outprefix

        if self.outdir is None:
            self.outdir = join(TMPDIR, 'tmp.mustache.minimus.'+str(randint(0, 1e20)))
        if self.outprefix is None:
            self.outprefix = 'mustache.minimus.' + str(randint(0, 1e20)) +'.' + \
                             str(hash(str(sample(reads, len(reads)))))

        self.reads = reads
        self.quals = quals
        self.read_names = read_names

        self.full_outprefix = join(self.outdir, self.outprefix)
        self.reads_path = self.full_outprefix + '.reads.fasta'
        self.quals_path = self.full_outprefix + '.quals.fasta'
        self.align_seq_fasta_path = self.full_outprefix + '.align_seq.fasta'
        self.align_sam_path = self.full_outprefix + '.align_seq.sam'
        self.afg_path = self.full_outprefix + '.afg'
        self.bam_path = self.full_outprefix + '.bam'

        os.makedirs(self.outdir, exist_ok=True)


    def delete_files(self):
        shell('rm -rf {reads} {quals} {afg} {align_seq_fasta} {align_sam} '
              '{outprefix}.bnk {outprefix}.fasta'.format(outprefix=self.full_outprefix, afg=self.afg_path,
                                                         align_seq_fasta=self.align_seq_fasta_path,
                                                         align_sam=self.align_sam_path, quals=self.quals_path,
                                                         reads=self.reads_path))


    def delete_afg_bank(self):
        shell('rm -rf {afg} {outprefix}.bnk'.format(outprefix=self.full_outprefix, afg=self.afg_path))


    def assemble(self, min_overlap_length=10, min_overlap_count=2, stranded=False):
        self.write_reads_as_fasta()
        self.delete_afg_bank()

        #print("toAmos")
        bank_file = self.full_outprefix + '.bnk'
        if self.quals:
            while not (isfile(self.reads_path) and isfile(self.quals_path)):
                continue

            shell("TMPDIR={tmpdir}; toAmos -s \"{reads_path}\" -q {quals_path} -o \"{afg_path}\" &> /dev/null;".format(
                tmpdir=self.outdir, reads_path=self.reads_path, quals_path=self.quals_path, afg_path=self.afg_path), read=True)
        else:
            while not (isfile(self.reads_path) and isfile(self.quals_path)):
                continue

            shell("toAmos -s \"{reads_path}\" -o \"{afg_path}\" &> /dev/null;".format(
                reads_path=self.reads_path, afg_path=self.afg_path), read=True)
        #print("bank-transact")
        while not isfile(self.afg_path):
            continue

        shell("bank-transact -f -z -b \"{bank_file}\" -m \"{afg_path}\" &> /dev/null;".format(
            bank_file=bank_file, afg_path=self.afg_path), read=True)
        #print("hash-overlap")

        if stranded:
            shell(
                "hash-overlap -s -o {min_olen} -B \"{bank_file}\" &> /dev/null;".format(bank_file=bank_file,
                                                                                        min_olen=min_overlap_length),
                read=True)
        else:
            shell("hash-overlap -o {min_olen} -B \"{bank_file}\" &> /dev/null;".format(
                bank_file=bank_file, min_olen=min_overlap_length), read=True)

        #print("tigger")
        shell("tigger -b \"{bank_file}\" &> /dev/null;".format(bank_file=bank_file), read=True)

        #print("make-consensus")
        shell("make-consensus -o {min_ocount} -B -b \"{bank_file}\" &> /dev/null".format(
            bank_file=bank_file, min_ocount=min_overlap_count), read=True)

        #print("bank2fasta")
        self.out_fasta = join(self.full_outprefix + '.fasta')
        shell("bank2fasta -d -b \"{bank_file}\" > \"{out_fasta}\" 2> /dev/null;".format(
            bank_file=bank_file, out_fasta=self.out_fasta), read=True)

        while not isfile(self.out_fasta):
            continue

        return self.out_fasta


    def write_reads_as_fasta(self):
        with open(self.reads_path, 'w') as out:
            for i in range(len(self.reads)):
                if self.read_names:
                    out.write('>%s' % self.read_names[i] + '\n' + self.reads[i] + '\n')
                else:
                    out.write('>Read%d' % i + '\n' + self.reads[i] + '\n')

        if self.quals:
            with open(self.quals_path, 'w') as out:
                for i in range(len(self.quals)):
                    outquals = ' '.join(map(str, query_qualities_to_phred(self.quals[i])))
                    if self.read_names:
                        out.write('>%s' % self.read_names[i] + '\n' + outquals + '\n')
                    else:
                        out.write('>Read%d' % i + '\n' + outquals + '\n')


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


    def something_assembled(self):
        for contig in SeqIO.parse(self.out_fasta, 'fasta'):
            return True
        return False


    def write_seq_to_fasta(self, seq):
        with open(self.align_seq_fasta_path, 'w') as out:
           out.write('>align_seq\n' + seq)


    def align_seq_to_assembly(self, seq):
        self.write_seq_to_fasta(seq)
        index_genome(self.out_fasta)
        align_to_genome_se(self.align_seq_fasta_path, self.out_fasta, self.align_sam_path)
        shell('rm -f %s.*' % self.out_fasta)


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
                    current_contig = str(contig.seq)

            if orient == 'R' and (not read.is_reverse):
                if len(current_contig) - read.reference_start <= query_length:
                    return None
                softclip = left_softclipped_sequence_strict(read)
                extension = softclip + current_contig[read.reference_start:]
            elif orient == 'R' and read.is_reverse:
                if read.reference_end <= query_length:
                    return None
                softclip = right_softclipped_sequence_strict(read)
                extension = revcomp(current_contig[:read.reference_end] + softclip)
            elif orient == 'L' and (not read.is_reverse):
                if read.reference_end <= query_length:
                    return None
                softclip = right_softclipped_sequence_strict(read)
                extension = current_contig[:read.reference_end] + softclip
            elif orient == 'L' and read.is_reverse:
                if len(current_contig) - read.reference_start <= query_length:
                    return None
                softclip = left_softclipped_sequence_strict(read)
                extension = revcomp(softclip + current_contig[read.reference_start:])

        return extension
