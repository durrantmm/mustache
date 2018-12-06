import warnings
warnings.filterwarnings("ignore")
import sys
import click
import pysam
import pygogo as gogo
import pandas as pd
import numpy as np
from snakemake import shell
from random import randint
from mustache import fastatools, embosstools, sctools, pysamtools, minimustools
from mustache.misc import revcomp
from os.path import basename
from multiprocessing import Pool
from os.path import join, dirname

verbose=True
logger = gogo.Gogo(__name__, verbose=verbose).logger

def extend(row):
    bam = pysam.AlignmentFile(row['bam_path'], 'rb')
    contig, pos_5p, pos_3p, seq_5p, seq_3p = row['contig'], row['pos_5p'], row['pos_3p'], row['seq_5p'], row['seq_3p']
    tmp_outdir = row['outdir']

    logger.info("Attempting to extend flank sequences at %s..." % ' '.join(map(str, [contig, pos_3p, pos_5p])))

    extended_seq_5p = get_extended_sequence(bam, contig, pos_5p, seq_5p, 'R', tmp_outdir)
    extended_seq_3p = get_extended_sequence(bam, contig, pos_3p, seq_3p, 'L', tmp_outdir)

    if extended_seq_5p != seq_5p:
        logger.info("5' sequence at %s extended successfully by %d base pairs" % (
            ':'.join(map(str, [contig, pos_5p])), len(extended_seq_5p) - len(seq_5p)))
    if extended_seq_3p != seq_3p:
        logger.info("3' sequence at %s extended successfully by %d base pairs" % (
            ':'.join(map(str, [contig, pos_3p])), len(extended_seq_3p) - len(seq_3p)))

    return extended_seq_5p, extended_seq_3p

def get_extended_sequence(bam, contig, pos, seq, orient, tmp_outdir):
    seq = seq.upper()
    reads, quals = get_reads_to_assemble(bam, contig, pos, orient, get_quals=True)

    if len(reads) == 0:
        return seq

    # print("RUNNING ASSEMBLY")
    assembler = minimustools.MinimusAssembler(reads, quals, outdir=tmp_outdir)
    assembler.assemble()

    # print("CHECKING IF SOMETHING ASSEMBLED")
    if not assembler.something_assembled():
        assembler.delete_files()
        return seq

    # print("ALIGNING TO ASSEMBLY")
    assembler.align_seq_to_assembly(seq)
    # print("RETRIEVING EXTENDED SEQUENCE")
    extended_seq = assembler.retrieve_extended_sequence(orient)
    assembler.delete_files()

    if extended_seq is None:
        return seq

    return extended_seq.upper()


def get_reads_to_assemble(bam, contig, pos, orient, get_quals=False):

    if get_quals:

        if orient == 'R':
            softclipped_reads, softclipped_quals = pysamtools.get_right_softclipped_reads_at_site(bam, contig, pos, get_quals=True)
            unmapped_reads, unmapped_quals = pysamtools.get_right_unmapped_reads(bam, contig, pos, get_quals=True)
        elif orient == 'L':
            softclipped_reads, softclipped_quals = pysamtools.get_left_softclipped_reads_at_site(bam, contig, pos, get_quals=True)
            unmapped_reads, unmapped_quals = pysamtools.get_left_unmapped_reads(bam, contig, pos, get_quals=True)

    else:
        if orient == 'R':
            softclipped_reads = pysamtools.get_right_softclipped_reads_at_site(bam, contig, pos, get_quals=True)
            unmapped_reads = pysamtools.get_right_unmapped_reads(bam, contig, pos, get_quals=True)
        elif orient == 'L':
            softclipped_reads = pysamtools.get_left_softclipped_reads_at_site(bam, contig, pos, get_quals=True)
            unmapped_reads = pysamtools.get_left_unmapped_reads(bam, contig, pos, get_quals=True)


    if get_quals:
        return softclipped_reads + unmapped_reads, softclipped_quals + unmapped_quals
    else:
        return softclipped_reads + unmapped_reads


def _extendpairs(pairsfile, bamfile, threads, output_file):
    pairs = pd.read_csv(pairsfile, sep='\t')

    if pairs.shape[0] == 0:
        logger.info("No pairs found in the input file...")

    else:
        sequences_5p = list(pairs['seq_5p'])
        sequences_3p = list(pairs['seq_3p'])

        logger.info("Running extendpairs algorithm on %d total pairs..." % pairs.shape[0])

        pair_rows = [row for index, row in pairs.iterrows()]
        tmp_outdir = join(dirname(output_file), 'tmp.mustache.minimus.' + str(randint(1, 1e20)))
        for row in pair_rows:
            row['bam_path'] = bamfile
            row['outdir'] = tmp_outdir

        agents = threads
        with Pool(processes=agents) as pool:
            extensions = np.array(pool.map(extend, pair_rows))
        shell("rmdir %s" % tmp_outdir)

        extensions_5p = [pair[0] for pair in extensions]
        extensions_3p = [pair[1] for pair in extensions]

        pairs['extended_5p'] = list(pairs['seq_5p'] != extensions_5p)
        pairs['extended_3p'] = list(pairs['seq_3p'] != extensions_3p)
        pairs['seq_5p'] = extensions_5p
        pairs['seq_3p'] = extensions_3p

    if output_file:
        logger.info("Saving results to file %s" % output_file)
        pairs.to_csv(output_file, sep='\t', index=False)

    return pairs


if __name__ == '__main__':
    extendpairs()