import sys
from mustache import output
from collections import OrderedDict
import pandas as pd
from os.path import join
from mustache.bam_tools import *
from mustache.bwa_tools import *
from mustache.assemble_tools import MinimusAssembler, Minimus2Merger

pd.set_option('display.max_rows', 1000)
pd.set_option('display.max_columns', 1000)
pd.set_option('display.width', 1000)


def process_candidate_site(bam_file, genome, contig, site, orientation, flank_length, max_mapping_distance,
                           average_read_length, assembly_kmer_size, outdir, output_prefix):

    site_name = contig + ':' + str(site)+ '_' + orientation

    click.echo('Processing the site %s' % site_name)
    bam_path = bam_file.filename.decode('UTF-8')
    bam_file = pysam.AlignmentFile(bam_path, 'rb')

    # Get the relevant reads
    softclipped_reads, unmapped_reads = None, None
    if orientation == "L":
        softclipped_reads = get_left_softclipped_reads_at_site(bam_file, contig, site, revcomp_left=True)
        unmapped_reads = get_left_unmapped_reads_in_flanks(bam_file, contig, site, flank_length, max_mapping_distance)
    elif orientation == "R":
        softclipped_reads = get_right_softclipped_reads_at_site(bam_file, contig, site, revcomp_left=True)
        unmapped_reads = get_right_unmapped_reads_in_flanks(bam_file, contig, site, flank_length, max_mapping_distance)

    # Get Flank Assembly
    flank_assembly = get_flank_assembly(bam_file, contig, site, orientation, softclipped_reads,
                                        unmapped_reads, outdir, site_name)
    if not flank_assembly:
        click.echo("\tFlank assembly failed...")

        return None

    click.echo("\tFlank Assembly:")
    click.echo(misc.wrap_string(flank_assembly, newline_char='\t\t'))

    # Get local ancestral genome and realign reads
    click.echo('\tRealigning reads to reconstructed ancestral genome and counting alignments...')
    junction_read_count, total_read_count = get_flank_read_count(genome, contig, site, orientation, flank_assembly, unmapped_reads,
                                            softclipped_reads, flank_length, outdir, site_name)

    click.echo('\tTotal reads aligning to assembled flank: %d' % total_read_count)

    out_dict = OrderedDict([
        ('contig', [contig]),
        ('site', [site]),
        ('orientation', [orientation]),
        ('flank_length', [flank_length]),
        ('inseq_read_count', [total_read_count]),
        ('assembly_length', [len(flank_assembly)]),
        ('flank_assembly', [flank_assembly])
    ])

    return pd.DataFrame(out_dict)


def get_flank_assembly(bam_file, contig, site, orientation, softclipped_reads, unmapped_reads, outdir, site_name):

    click.echo('\tPerforming assembly...')

    tmp_sc_fasta = join(outdir, site_name + '.sc.flank_assembly.fasta')
    tmp_sc_alignment = join(outdir, site_name + ".sc.flank_assembly.bam")

    click.echo('\tRunning minimus assembler...')
    assembler = MinimusAssembler([seq[1] for seq in softclipped_reads + unmapped_reads], outdir, site_name)
    output_assembly = assembler.assemble()

    contig_count = 0
    for rec in SeqIO.parse(output_assembly, 'fasta'):
        contig_count += 1
    if contig_count == 0:
        return None
    click.echo('\t%d assembled contigs available at %s' % (contig_count, output_assembly))

    click.echo('\tAligning softclipped ends to find matching contig...')
    output.write_reads_to_fasta(softclipped_reads, tmp_sc_fasta, seq_index=2)

    index_genome(output_assembly, silence=True)
    bwa_aln_to_genome_single(tmp_sc_fasta, output_assembly, tmp_sc_alignment, silence=True)

    click.echo('\tFinding exact site of insertion in contig...')
    flank_assembly = retrieve_flanking_sequence(output_assembly, pysam.AlignmentFile(tmp_sc_alignment, 'r'), orientation)

    click.echo('\tDeleting tmp files...')
    shell("rm -f {output_assembly}* {tmp_sc_fasta} {tmp_sc_alignment}*;")
    assembler.delete_files()

    return flank_assembly


def determine_site_pairs(df, bam_file, genome_fasta, min_pair_distance, max_pair_distance, outdir):

    click.echo("\nDetermining insertion site pairs and attempting to merge flanks...")

    pairwise_info = None
    genome_name = None
    genome_seq = None

    for index1, left_mate in df.iterrows():

        if left_mate['orientation'] == 'R':
            continue

        contig = left_mate['contig']
        left_site = left_mate['site']
        left_inseq_count = left_mate['inseq_read_count']
        left_flank = left_mate['flank_assembly']

        # Retrieve the relevant contig from the fasta
        if genome_name != contig:
            click.echo("\tLoading relevant genome from fasta...")
            read_genome_fasta = SeqIO.parse(genome_fasta, format='fasta')

            for rec in read_genome_fasta:
                if rec.name == contig:
                    genome_name = rec.name
                    genome_seq = str(rec.seq)
                    break

        nearby_mates = misc.get_nearby_sites(df, left_mate, min_pair_distance, max_pair_distance)
        for index2, right_mate in nearby_mates.iterrows():

            right_site = right_mate['site']
            right_inseq_count = right_mate['inseq_read_count']
            right_flank = right_mate['flank_assembly']

            site_name = contig + ':' + str(left_site) + '-' + str(right_site)

            total_count = left_inseq_count + right_inseq_count
            count_diff = abs(left_inseq_count - right_inseq_count)

            left_softclipped_reads = get_left_softclipped_reads_at_site(bam_file, contig, left_site, revcomp_left=True)
            right_softclipped_reads = get_right_softclipped_reads_at_site(bam_file, contig, right_site, revcomp_left=True)

            merged_assembly = attempt_flank_merge(genome_seq, contig, left_site, right_site, left_flank, right_flank,
                                                  left_softclipped_reads, right_softclipped_reads, outdir, site_name)

            merged_length = -1
            if merged_assembly:
                merged_length = len(merged_assembly)

            pair_dict = OrderedDict([
                ('contig', [contig]),
                ('left_site', [int(left_site)]),
                ('right_site', [int(right_site)]),
                ('left_inseq_read_count', [left_inseq_count]),
                ('right_inseq_read_count', [right_inseq_count]),
                ('total_count', [total_count]),
                ('count_diff', [count_diff]),
                ('merged_assembly_length', [merged_length]),
                ('merged_assembly', [merged_assembly]),
            ])

            pair_df = pd.DataFrame(pair_dict)
            if pairwise_info is None:
                pairwise_info = pair_df
            else:
                pairwise_info = pd.concat([pairwise_info, pair_df ])

    return pairwise_info


def attempt_flank_merge(genome, contig, left_site, right_site, left_flank, right_flank,
                        left_softclipped_reads, right_softclipped_reads, outdir, site_name, assembly_kmer_size=11):

    tmp_left_sc_fasta = join(outdir, site_name + ".sc.left.flank_merge.fasta")
    tmp_right_sc_fasta = join(outdir, site_name + ".sc.right.flank_merge.fasta")
    tmp_left_sc_alignment = join(outdir, site_name + ".sc.left.flank_merge.bam")
    tmp_right_sc_alignment = join(outdir, site_name + ".sc.right.flank_merge.bam")

    left_flank = get_inserted_genome(genome, contig, left_site, 'L', left_flank, 200)
    right_flank = get_inserted_genome(genome, contig, right_site, 'R', right_flank, 200)

    output.write_reads_to_fasta(left_softclipped_reads, tmp_left_sc_fasta, seq_index=2)
    output.write_reads_to_fasta(right_softclipped_reads, tmp_right_sc_fasta, seq_index=2)

    merger = Minimus2Merger([left_flank, right_flank], outdir, site_name)
    output_assembly = merger.merge(return_assembly_path=True)

    if not output_assembly:
        shell("rm -f {output_assembly}* {tmp_left_sc_fasta} {tmp_right_sc_fasta} "
              "{tmp_left_sc_alignment}* {tmp_right_sc_alignment}*;")
        merger.delete_files()
        return None

    index_genome(output_assembly, silence=True)
    bwa_aln_to_genome_single(tmp_left_sc_fasta, output_assembly, tmp_left_sc_alignment, silence=True)
    bwa_aln_to_genome_single(tmp_right_sc_fasta, output_assembly, tmp_right_sc_alignment, silence=True)

    merged_assembly = retrieve_merged_assembly(output_assembly,
                                               pysam.AlignmentFile(tmp_left_sc_alignment, 'rb'),
                                               pysam.AlignmentFile(tmp_right_sc_alignment, 'rb'))

    shell("rm -f {output_assembly}* {tmp_left_sc_fasta} {tmp_right_sc_fasta} "
          "{tmp_left_sc_alignment}* {tmp_right_sc_alignment}*;")
    merger.delete_files()

    return merged_assembly


def prioritize_site_pairs(site_pairs):
    keep_pairs = set()
    keep_sites = set()

    site_pairs = site_pairs.sort_values(['merged_assembly_length', 'count_diff', 'total_count'], ascending=False).reset_index(drop=True)

    for index, row in site_pairs.iterrows():
        contig, left_site, right_site = row['contig'], row['left_site'], row['right_site']

        if (contig, left_site) in keep_sites or (contig, right_site) in keep_sites:
            continue

        keep_sites.add((contig, left_site))
        keep_sites.add((contig, right_site))

        keep_pairs.add((contig, left_site, right_site))

    final_keep_pairs = None
    for contig, left_site, right_site in keep_pairs:
        pair = site_pairs.query('contig == "{contig}" & left_site == {left_site} & right_site == {right_site}'.format(
            contig=contig, left_site=left_site, right_site=right_site
        ))

        if final_keep_pairs is None:
            final_keep_pairs = pair
        else:
            final_keep_pairs = pd.concat([final_keep_pairs, pair])

    return final_keep_pairs


def create_final_output_table(assembled_flanks, priority_site_pairs):

    final_table = None
    for index, row in priority_site_pairs.iterrows():
        contig, left_site, right_site = row['contig'], int(row['left_site']), int(row['right_site'])
        left_inseq_count, right_inseq_count = int(row['left_inseq_read_count']), int(row['right_inseq_read_count'])
        merged_assembly = row['merged_assembly']

        if merged_assembly:
            out_dict = OrderedDict([
                ('contig', [contig]),
                ('site', [left_site]),
                ('left_site', [left_site]),
                ('right_site', [right_site]),
                ('orientation', ['M']),
                ('partner_site', [-1]),
                ('left_site_read_count', [left_inseq_count]),
                ('right_site_read_count', [right_inseq_count]),
                ('assembly_length', [len(merged_assembly)]),
                ('assembly', [merged_assembly])
            ])
        else:
            left_flank_assembly = assembled_flanks.query("site == @left_site").reset_index(drop=True).at[0,'flank_assembly']
            right_flank_assembly = assembled_flanks.query("site == @right_site  ").reset_index(drop=True).at[0,'flank_assembly']

            out_dict = OrderedDict([
                ('contig', [contig, contig]),
                ('site', [left_site, right_site]),
                ('left_site', [left_site, -1]),
                ('right_site', [-1, right_site]),
                ('orientation', ['L','R']),
                ('partner_site', [right_site, left_site]),
                ('left_site_read_count', [left_inseq_count, -1]),
                ('right_site_read_count', [-1, right_inseq_count]),
                ('assembly_length', [len(left_flank_assembly), len(right_flank_assembly)]),
                ('assembly', [left_flank_assembly, right_flank_assembly])
            ])

        if final_table is None:
            final_table = pd.DataFrame(out_dict)
        else:
            final_table = pd.concat([final_table, pd.DataFrame(out_dict)])

    return final_table.sort_values(['contig', 'site']).reset_index(drop=True).drop(columns=['site'])

