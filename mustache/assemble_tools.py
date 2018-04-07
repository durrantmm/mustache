import sys
from snakemake import shell
import os
from os.path import join


class MinimusAssembler:

    full_outprefix = None

    reads_path = None
    reads = None
    afg_path = None
    out_fasta = None

    def __init__(self, reads, outdir='tmp', outprefix='tmp'):

        self.reads = reads

        self.full_outprefix = join(outdir, outprefix)
        self.reads_path = join(outdir, outprefix + '.reads.fasta')
        self.afg_path = join(outdir, outprefix + '.afg')

        os.makedirs(outdir, exist_ok=True)

    def delete_files(self):
        shell('rm -rf {reads} {afg} {outprefix}.bnk {outprefix}.fasta'.format(outprefix=self.full_outprefix,
                                                                              afg=self.afg_path, reads=self.reads_path))

    def assemble(self, min_overlap_length=10, min_overlap_count=2):
        self.write_reads_as_fasta()
        shell("toAmos -s {reads_path} -o {afg_path} 2> junk;".format(reads_path=self.reads_path, afg_path=self.afg_path), read=True)
        shell("bank-transact -f -z -b {out_prefix}.bnk -m {out_prefix}.afg 2> junk;".format(out_prefix=self.full_outprefix), read=True)
        shell("hash-overlap -o {min_olen} -B {out_prefix}.bnk 2> junk;".format(out_prefix=self.full_outprefix, min_olen=min_overlap_length), read=True)
        shell("tigger -b {out_prefix}.bnk 2> junk;".format(out_prefix=self.full_outprefix), read=True)
        shell("make-consensus -o {min_ocount} -B -b {out_prefix}.bnk 2> junk;".format(out_prefix=self.full_outprefix, min_ocount=min_overlap_count), read=True)
        shell("bank2fasta -b {out_prefix}.bnk > {out_prefix}.fasta 2> junk;".format(out_prefix=self.full_outprefix), read=True)
        shell('rm -f junk;')

        self.out_fasta = join(self.full_outprefix + '.fasta')
        return self.out_fasta

    def write_reads_as_fasta(self):
        with open(self.reads_path, 'w') as out:
            for i in range(len(self.reads)):
                out.write('>Read%d' % i + '\n' + self.reads[i] + '\n')



class Minimus2Merger:

    full_outprefix = None

    reads_path = None
    reads = None
    afg_path = None
    out_fasta = None

    def __init__(self, reads, outdir='tmp', outprefix='tmp'):

        self.reads = reads

        self.full_outprefix = join(outdir, outprefix)
        self.reads_path = join(outdir, outprefix + '.reads.fasta')
        self.afg_path = join(outdir, outprefix + '.afg')

        os.makedirs(outdir, exist_ok=True)

    def delete_files(self):
        shell('rm -rf {reads} {afg} {outprefix}.bnk {outprefix}.fasta'.format(outprefix=self.full_outprefix,
                                                                              afg=self.afg_path, reads=self.reads_path))

    def assemble(self):
        self.write_reads_as_fasta()

        shell("toAmos -s {reads_path} -o {afg_path} 2> junk;".format(reads_path=self.reads_path, afg_path=self.afg_path),read=True)
        shell("bank-transact -c -z -b {outprefix}.bnk -m {afg_path};".format(outprefix=self.full_outprefix, afg_path=self.afg_path), read=True)
        shell("dumpreads {outprefix}.bnk -M 1 > {outprefix}.ref.seq;".format(outprefix=self.full_outprefix), read=True)
        shell("dumpreads {outprefix}.bnk -m 1 > {outprefix}.qry.seq;".format(outprefix=self.full_outprefix), read=True)
        shell("nucmer -maxmatch -c 10 {outprefix}.ref.seq {outprefix}.qry.seq -p {outprefix};".format(outprefix=self.full_outprefix), read=True)
        shell("show-coords -H -c -l -o -r -I 94 {outprefix}.delta | nucmerAnnotate | egrep 'BEGIN|END|CONTAIN|IDENTITY' > {outprefix}.coords;".format(outprefix=self.full_outprefix))
        shell("nucmer2ovl -ignore 20 -tab {outprefix}.coords | sort2 > {outprefix}.ovl;".format(outprefix=self.full_outprefix), read=True)
        shell("ovl2OVL {outprefix}.ovl  > {outprefix}.OVL;".format(outprefix=self.full_outprefix), read=True)
        shell("bank-transact -z -b {outprefix}.bnk -m {outprefix}.OVL;".format(outprefix=self.full_outprefix), read=True)
        shell("tigger -b {outprefix}.bnk;".format(outprefix=self.full_outprefix), read=True)
        shell("make-consensus -B -e 0.06 -b {outprefix}.bnk -w 15;".format(outprefix=self.full_outprefix), read=True)
        shell("bank2fasta -b {outprefix}.bnk > {outprefix}.fasta;".format(outprefix=self.full_outprefix), read=True)

        self.out_fasta = join(self.full_outprefix + '.fasta')
        return self.out_fasta

    def write_reads_as_fasta(self):
        with open(self.reads_path, 'w') as out:
            for i in range(len(self.reads)):
                out.write('>Read%d' % i + '\n' + self.reads[i] + '\n')