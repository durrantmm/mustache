import sys
from os.path import join
import click
from mustache import output, bwa_tools, misc, bam_tools, bowtie2
import pysam
from collections import defaultdict, OrderedDict
from Bio import SeqIO
import numpy as np
import pandas as pd
from snakemake import shell

def run(insertion_seqs, comp_genome, output_prefix, outdir):

    click.echo("Creating insertion sequence record names")
    insertion_seqs = create_record_names(insertion_seqs)

    num_merged = insertion_seqs.query('orientation=="M"').shape[0]
    num_flanks = insertion_seqs.query('orientation!="M"').shape[0]

    merged_outfasta = output.write_dataframe_to_fasta_by_columns(insertion_seqs.query('orientation=="M"'), 'name', 'assembly', outdir, 'merged')
    flanks_outfasta1, flanks_outfasta2 = output.write_flanks_to_paired_end_fasta_by_columns(insertion_seqs.query('orientation!="M"'), 'name', 'assembly', outdir, 'flanks')

    click.echo("Indexing genome if not already indexed...")
    is_indexed = bowtie2.genome_is_indexed(comp_genome)
    if not is_indexed:
        bowtie2.index_genome(comp_genome, silence=True)
    click.echo("Genome has been indexed...")


    merged_bam = join(outdir, 'merged.realigned.bam')
    if num_merged > 0:
        bowtie2.align_fasta_to_genome(merged_outfasta, comp_genome, merged_bam, threads=1, silence=True)

    flanks_bam = join(outdir, 'flanks.realigned.bam')
    if num_flanks > 0:
        bowtie2.align_paired_fasta_to_genome(flanks_outfasta1, flanks_outfasta2, comp_genome, flanks_bam, threads=1,
                                             silence=True, additional_flags='-X 1000000 --no-mixed --no-discordant')
    #sys.exit()
    click.echo("Recovering the mapped reference for full sequences...")
    if num_merged > 0:
        insertion_seqs = recover_merged_sequences(insertion_seqs, merged_bam, comp_genome)
    else:
        insertion_seqs['full_recovery'] = np.nan

    click.echo("Recovering the mapped reference for flank sequences...")
    if num_flanks > 0:
        insertion_seqs = recover_merged_flank_sequences(insertion_seqs, flanks_bam, comp_genome)
    else:
        insertion_seqs['flank_recovery'] = np.nan

    insertion_seqs['assembly'] = insertion_seqs.apply(lambda r: r.recovered_assembly.upper() if r.flank_recovery is True or r.full_recovery is True else r.assembly.upper(), axis=1)
    insertion_seqs['assembly_length'] = insertion_seqs.apply(lambda r: len(r.assembly), axis=1)
    final_data = pd.DataFrame(insertion_seqs[['contig', 'left_site', 'right_site', 'orientation', 'partner_site',
                                              'full_recovery', 'flank_recovery', 'assembly_length', 'assembly']])

    click.echo("Writing final %d insertion sequence records to file..." % final_data.shape[0])
    output.write_final_dataframe(final_data, outdir, output_prefix, suffix='.recovered_insertion_seqs.tsv')
    output.write_final_merged_dataframe_to_fasta(final_data, outdir, output_prefix, suffix='.recovered_insertion_seqs.fasta')

    click.echo("Deleting temporary files...")
    shell("rm {file1}* {file2}* {file3}* {file4}* {file5}*".format(file1=merged_outfasta, file2=flanks_outfasta1,
                                                                   file3=flanks_outfasta2, file4=merged_bam,
                                                                   file5=flanks_bam))


def recover_merged_sequences(insertion_seqs, bam_file, genome_fasta):
    out_seqs = pd.DataFrame(insertion_seqs)
    out_seqs['full_recovery'] = np.nan

    bam = pysam.AlignmentFile(bam_file, 'rb')

    for read in bam:
        record_index = int(read.query_name.split('|')[-1])

        if bam_tools.is_left_softclipped_strict(read) or bam_tools.is_right_softclipped_strict(read):
            continue

        recovered_assembly = read.get_reference_sequence().upper()

        if read.is_reverse:
            recovered_assembly = misc.revcomp(recovered_assembly)


        out_seqs.loc[record_index, 'recovered_assembly'] = recovered_assembly
        out_seqs.loc[record_index, 'full_recovery'] = True

    return out_seqs


def recover_merged_flank_sequences(insertion_seqs, bam_file, genome_fasta):
    out_seqs = pd.DataFrame(insertion_seqs)
    out_seqs['flank_recovery'] = np.nan

    bam = pysam.AlignmentFile(bam_file, 'rb')

    contigs = misc.fastaContigGetter(genome_fasta)

    merged_flank_indices = []

    for read1 in bam:
        genome = contigs.get_contig(read1.reference_name)

        if read1.is_reverse:
           continue

        read2 = bam.mate(read1)


        if read1.is_read1 and read2.is_read2:
            right_read, left_read = read1, read2
        else:
            right_read, left_read = read2, read1

        left_index = int(left_read.query_name.split('||')[0].split('|')[-1])
        right_index = int(right_read.query_name.split('||')[-1].split('|')[-1])

        merged_flank_indices.append((left_index, right_index))

        if not right_read.is_reverse:

            right_read_softclip = bam_tools.get_left_softclip_length_strict(right_read)
            left_read_softclip = bam_tools.get_right_softclip_length_strict(left_read)

            right_read_softclipping = right_read.query_sequence[:right_read_softclip]
            left_read_softclipping = left_read.query_sequence[-left_read_softclip:]

            if right_read_softclip == 0:
                right_read_softclipping = ''
            if left_read_softclip == 0:
                left_read_softclipping = ''

            #if read1.query_name == 'NC_000913.3|69096|-1|L|2||NC_000913.3|-1|69096|R|3':
            #   print(right_read.cigartuples)
            #    print(left_read.cigartuples)
            #    print(right_read_softclip)
            #    print(left_read_softclip)
            #    sys.exit()

            right_read_start = right_read.reference_start
            left_read_end = left_read.reference_end

            recovered_assembly = right_read_softclipping + genome[right_read_start:left_read_end] + left_read_softclipping

            out_seqs.loc[left_index, 'recovered_assembly'] = recovered_assembly
            out_seqs.loc[right_index, 'recovered_assembly'] = recovered_assembly
            out_seqs.loc[left_index, 'flank_recovery'] = True
            out_seqs.loc[right_index, 'flank_recovery'] = True

        else:


            right_read_softclip = bam_tools.get_right_softclip_length_strict(right_read)
            left_read_softclip = bam_tools.get_left_softclip_length_strict(left_read)

            right_read_softclipping = right_read.query_sequence[-right_read_softclip:]
            left_read_softclipping = left_read.query_sequence[:left_read_softclip]

            if right_read_softclip == 0:
                right_read_softclipping = ''
            if left_read_softclip == 0:
                left_read_softclipping = ''


            right_read_end = right_read.reference_end
            left_read_start = left_read.reference_start

            recovered_assembly = left_read_softclipping + genome[left_read_start:right_read_end] + right_read_softclipping

            recovered_assembly = misc.revcomp(recovered_assembly.upper())

            out_seqs.loc[left_index, 'recovered_assembly'] = recovered_assembly
            out_seqs.loc[right_index, 'recovered_assembly'] = recovered_assembly
            out_seqs.loc[left_index, 'flank_recovery'] = True
            out_seqs.loc[right_index, 'flank_recovery'] = True

    for i1, i2 in merged_flank_indices:
        out_seqs.loc[i1, 'right_site'] = out_seqs.loc[i2, 'right_site']
        out_seqs.loc[i1, 'orientation'] = 'M'
        out_seqs.loc[i1, 'partner_site'] = -1

    out_seqs = out_seqs.query('orientation != "R" | flank_recovery != True')
    return out_seqs


def parse_read_name(read_name):
    if read_name.find("MERGED") != -1:
        contig, left_site, right_site, orientation = read_name.split('|')[:-2]

    elif read_name.find("LEFT") != -1:
        contig, left_site, orientation, partner_name, right_site = read_name.split('|')[:-2]

    else:
        contig, right_site, orientation, partner_name, left_site = read_name.split('|')[:-2]

    left_site, right_site = int(left_site), int(right_site)

    return contig, left_site, right_site, orientation


def create_record_names(insertion_seqs):
    out_seqs = pd.DataFrame(insertion_seqs)
    out_names = []
    for index, row in insertion_seqs.iterrows():
        out_names.append('|'.join(map(str, [row.contig, row.left_site, row.right_site, row.orientation, index])))

    out_seqs['name'] = out_names

    return out_seqs