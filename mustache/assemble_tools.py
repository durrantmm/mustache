import sys
from snakemake import shell
import os
from os.path import join
import subprocess

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
    def delete_afg_bank(self):
        shell('rm -rf {afg} {outprefix}.bnk'.format(outprefix=self.full_outprefix, afg=self.afg_path, reads=self.reads_path))

    def assemble(self, min_overlap_length=10, min_overlap_count=2):
        self.write_reads_as_fasta()
        self.delete_afg_bank()
        shell("toAmos -s {reads_path} -o {afg_path} 2> /dev/null;".format(reads_path=self.reads_path, afg_path=self.afg_path), read=True)
        shell("bank-transact -f -z -b {out_prefix}.bnk -m {out_prefix}.afg 2> /dev/null;".format(out_prefix=self.full_outprefix), read=True)
        shell("hash-overlap -o {min_olen} -B {out_prefix}.bnk 2> /dev/null;".format(out_prefix=self.full_outprefix, min_olen=min_overlap_length), read=True)
        shell("tigger -b {out_prefix}.bnk 2> /dev/null;".format(out_prefix=self.full_outprefix), read=True)
        shell("make-consensus -o {min_ocount} -B -b {out_prefix}.bnk 2> /dev/null;".format(out_prefix=self.full_outprefix, min_ocount=min_overlap_count), read=True)
        shell("bank2fasta -b {out_prefix}.bnk > {out_prefix}.fasta 2> /dev/null;".format(out_prefix=self.full_outprefix), read=True)

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
        self.delete_files()

    def delete_files(self):
        shell('rm -rf {reads} {afg} {outprefix}.bnk {outprefix}.fasta '
              '{outprefix}.coords {outprefix}.ovl {outprefix}.OVL '
              '{outprefix}.qry.seq {outprefix}.ref.seq {outprefix}.delta'.format(outprefix=self.full_outprefix,
                                                                              afg=self.afg_path, reads=self.reads_path))
    def delete_int_files(self):
        shell('rm -rf {afg} {outprefix}.bnk {outprefix}.coords {outprefix}.ovl {outprefix}.OVL '
              '{outprefix}.qry.seq {outprefix}.ref.seq {outprefix}.delta'.format(outprefix=self.full_outprefix,
                                                                                 afg=self.afg_path, reads=self.reads_path))

    def merge(self, ignore_overhangs=2, min_overlap_length=10):
        self.write_reads_as_fasta()
        shell("toAmos -s {reads_path} -o {afg_path} 2> /dev/null;".format(reads_path=self.reads_path, afg_path=self.afg_path),read=True)
        shell("bank-transact -f -z -b {outprefix}.bnk -m {afg_path} 2> /dev/null;".format(outprefix=self.full_outprefix, afg_path=self.afg_path), read=True)
        shell("dumpreads {outprefix}.bnk -M 1 > {outprefix}.ref.seq 2> /dev/null;".format(outprefix=self.full_outprefix), read=True)
        shell("dumpreads {outprefix}.bnk -m 1 > {outprefix}.qry.seq 2> /dev/null;".format(outprefix=self.full_outprefix), read=True)
        shell("nucmer -maxmatch -c {min_olap} {outprefix}.ref.seq {outprefix}.qry.seq -p {outprefix} 2> /dev/null;".format(outprefix=self.full_outprefix, min_olap=min_overlap_length), read=True)

        # Will fail if there are no overlaps
        try:
            shell("show-coords -H -c -l -o -r -I 94 {outprefix}.delta | nucmerAnnotate | egrep 'BEGIN|END|CONTAIN|IDENTITY' > {outprefix}.coords 2> /dev/null;".format(outprefix=self.full_outprefix))
        except subprocess.CalledProcessError:
            return None

        shell("nucmer2ovl -ignore {ignore} -tab {outprefix}.coords | sort2 > {outprefix}.ovl 2> /dev/null;".format(outprefix=self.full_outprefix, ignore=ignore_overhangs), read=True)
        shell("ovl2OVL {outprefix}.ovl  > {outprefix}.OVL 2> /dev/null;".format(outprefix=self.full_outprefix), read=True)
        shell("bank-transact -z -b {outprefix}.bnk -m {outprefix}.OVL 2> /dev/null;".format(outprefix=self.full_outprefix), read=True)
        shell("tigger -b {outprefix}.bnk 2> /dev/null;".format(outprefix=self.full_outprefix), read=True)
        shell("make-consensus -B -e 0.06 -b {outprefix}.bnk -w 15 2> /dev/null;".format(outprefix=self.full_outprefix), read=True)
        shell("bank2fasta -b {outprefix}.bnk > {outprefix}.fasta 2> /dev/null;".format(outprefix=self.full_outprefix), read=True)

        self.out_fasta = join(self.full_outprefix + '.fasta')

        if os.stat(self.out_fasta).st_size == 0:
            return None
        else:
            return self.out_fasta

    def write_reads_as_fasta(self):
        with open(self.reads_path, 'w') as out:
            for i in range(len(self.reads)):
                out.write('>Read%d' % i + '\n' + self.reads[i] + '\n')