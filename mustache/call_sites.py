import sys
import click
from Bio import SeqIO
import pysam
from mustache.bam_tools import *
from collections import OrderedDict
import pandas as pd
from tqdm import tqdm

def call_sites(candidate_sites, bam_file, genome_fasta, max_mapping_distance, flank_length,
               min_alignment_overlap, min_alignment_overlap_count, outdir):
    click.echo("Calling insertion sequences at specified sites...")

    #candidate_sites.sort_values(['contig'], inplace=True)

    site_inseq_counts = get_all_inseq_read_counts(candidate_sites, bam_file, genome_fasta, max_mapping_distance,
                                                  flank_length, min_alignment_overlap, min_alignment_overlap_count,
                                                  outdir)
    site_ancestral_counts = get_all_ancestral_read_counts(candidate_sites, bam_file)
    print(site_inseq_counts)
    sys.exit()


def get_all_inseq_read_counts(candidate_sites, bam_file, genome_fasta, max_mapping_distance, flank_length,
                          min_alignment_overlap, min_alignment_overlap_count, outdir):

    orig_sites = candidate_sites[['contig', 'left_site', 'right_site', 'orientation', 'partner_site', 'assembly']].drop_duplicates()
    unique_sites = candidate_sites[['contig', 'left_site', 'right_site', 'orientation', 'assembly']].drop_duplicates()

    read_genome_fasta = SeqIO.parse(genome_fasta, format='fasta')
    genome_name = None

    site_inseq_counts = None

    with tqdm(total=unique_sites.shape[0], desc='GETTING EVIDENCE FOR INSERTION SEQUENCES') as bar:

        for index, row in unique_sites.iterrows():
            bar.update(1)
            contig = row['contig']
            orientation = row['orientation']
            assembly = row['assembly']

            # Retrieve the relevant contig from the fasta
            if genome_name != contig:
                next_rec = next(read_genome_fasta)
                while next_rec.name != contig:
                    next_rec = next(read_genome_fasta)
                genome_name = next_rec.name
                genome_seq = str(next_rec.seq)

            if orientation == 'L':
                site = row['left_site']
                df = get_inseq_read_count(bam_file, genome_seq, contig, site, orientation, assembly, max_mapping_distance,
                               flank_length, min_alignment_overlap, min_alignment_overlap_count, outdir)
            elif orientation == 'R':
                site = row['right_site']
                df = get_inseq_read_count(bam_file, genome_seq, contig, site, orientation, assembly, max_mapping_distance,
                               flank_length, min_alignment_overlap, min_alignment_overlap_count, outdir)
            else:
                left_site = row['left_site']
                right_site = row['right_site']
                df = pd.concat([
                    get_inseq_read_count(bam_file, genome_seq, contig, left_site, 'L', assembly, max_mapping_distance,
                              flank_length, min_alignment_overlap, min_alignment_overlap_count, outdir),
                    get_inseq_read_count(bam_file, genome_seq, contig, right_site, 'R', assembly, max_mapping_distance,
                              flank_length, min_alignment_overlap, min_alignment_overlap_count, outdir)
                ])

            if site_inseq_counts is None:
                site_inseq_counts = df
            else:
                site_inseq_counts = pd.concat([site_inseq_counts, df])

    site_inseq_counts.rename(columns={'flank_assembly': 'assembly'}, inplace=True)

    site_inseq_counts.loc[:, 'left_site'] = site_inseq_counts.apply(lambda r: r.site if r.orientation=='L' else -1, axis=1)
    site_inseq_counts.loc[:, 'right_site'] = site_inseq_counts.apply(lambda r: r.site if r.orientation == 'R' else -1, axis=1)

    site_inseq_counts = pd.DataFrame(site_inseq_counts[['contig', 'left_site', 'right_site', 'orientation', 'assembly', 'inseq_read_count']])
    left_inseq_counts = site_inseq_counts.query('orientation == "L" & right_site == -1')[['contig', 'left_site', 'assembly', 'inseq_read_count']]
    left_inseq_counts.rename(columns={'inseq_read_count': 'left_inseq_read_count'}, inplace=True)
    right_inseq_counts = site_inseq_counts.query('orientation == "R" & left_site == -1')[['contig', 'right_site', 'assembly', 'inseq_read_count']]
    right_inseq_counts.rename(columns={'inseq_read_count': 'right_inseq_read_count'}, inplace=True)

    orig_sites = pd.merge(orig_sites, left_inseq_counts, how='left')
    orig_sites = pd.merge(orig_sites, right_inseq_counts, how='left')
    orig_sites.fillna(-1, inplace=True)
    orig_sites.loc[:,'left_inseq_read_count'] = orig_sites.apply(lambda r: int(r.left_inseq_read_count), axis=1)
    orig_sites.loc[:, 'right_inseq_read_count'] = orig_sites.apply(lambda r: int(r.right_inseq_read_count), axis=1)

    print(orig_sites)
    sys.exit()
    return site_inseq_counts


def get_inseq_read_count(bam_file, genome, contig, site, orientation, assembly, max_mapping_distance,
                         flank_length, min_alignment_overlap, min_alignment_overlap_count, outdir):

    site_name = contig + ':' + str(site)+ '_' + orientation

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

    # Get local ancestral genome and realign reads
    inseq_read_count = get_flank_read_count(genome, contig, site, orientation, assembly, unmapped_reads,
                                            softclipped_reads, flank_length, outdir, site_name, min_alignment_overlap,
                                            min_alignment_overlap_count)

    out_dict = OrderedDict([
        ('contig', [contig]),
        ('site', [site]),
        ('orientation', [orientation]),
        ('flank_length', [flank_length]),
        ('inseq_read_count', [inseq_read_count]),
        ('assembly_length', [len(assembly)]),
        ('flank_assembly', [assembly])
    ])

    return pd.DataFrame(out_dict)

def get_all_ancestral_read_counts(candidate_sites, bam_file, max_softclip_distance = 20):

    for pu in bam_file.pileup(contig, right_site-1, right_site, truncate=True):
        for pr in pu.pileups:
            read = pr.alignment

            if read_runs_through_site(read, left_site, right_site):

                if is_right_softclipped(read):

                    site_clipped_read = clip_read_at_right_site(read, right_site)
                    align_score = alignment_tools.left_alignment_score(site_clipped_read, right_assembly)
                    if align_score > 0 and read.query_name not in keep_read_names:
                        continue

                if is_left_softclipped(read):

                    site_clipped_read = clip_read_at_left_site(read, left_site)
                    align_score = alignment_tools.right_alignment_score(site_clipped_read, left_assembly)

                    if align_score > 0 and read.query_name not in keep_read_names:
                        continue

                if read.query_name not in keep_read_names:
                    keep_read_names.add(read.query_name)
                    ancestral_reads.append(read)
