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
        shell('rm -rf {outprefix}*'.format(outprefix=self.full_outprefix))

    def assemble(self):
        self.write_reads_as_fasta()
        shell("toAmos -s {reads_path} -o {afg_path}".format(reads_path=self.reads_path, afg_path=self.afg_path), read=True)
        shell("minimus {out_prefix}".format(out_prefix=self.full_outprefix), read=True)

        self.out_fasta = join(self.full_outprefix + '.fasta')
        return self.out_fasta

    def write_reads_as_fasta(self):
        with open(self.reads_path, 'w') as out:
            for i in range(len(self.reads)):
                out.write('>Read%d' % i + '\n' + self.reads[i] + '\n')
