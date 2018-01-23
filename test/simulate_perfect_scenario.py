import click
import sys
from random import randint
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from snakemake import shell

@click.command()
@click.argument('output_prefix')
@click.option('--inseq_length', default=1500, help="The length of the insertion sequence")
@click.option('--genome_length', default=10000, help="The length of the insertion sequence.")
@click.option('--direct_repeat_length', default=5, help="The length of the direct repeat created by the insertion sequence.")
@click.option('--read_length', default=100, help="The length of each read pair")
@click.option('--fragment_length', default=500, help="The length of the fragment used.")
@click.option('--step_size', default=1, help="Specify the step size between fragments.")
def main(output_prefix, inseq_length, genome_length, direct_repeat_length, read_length, fragment_length, step_size):

    inseq = generate_random_sequence(inseq_length)
    ancestral_genome = generate_random_sequence(genome_length)

    inseq_genome = insert_sequence_center(ancestral_genome, inseq, direct_repeat_length)

    ancestral_reads = get_reads(ancestral_genome, read_length, fragment_length, step_size)
    inseq_reads = get_reads(inseq_genome, read_length, fragment_length, step_size)

    ancestral_genome_fasta_path = output_prefix + '.ancestral_genome.fna'
    write_sequence(ancestral_genome, 'ancestral_genome', ancestral_genome_fasta_path)

    inseq_genome_fasta_path = output_prefix + '.inseq_genome.fna'
    write_sequence(inseq_genome, 'inseq_genome', inseq_genome_fasta_path)

    ancestral_fastq_path = output_prefix + '.ancestral_reads'
    write_reads_to_fastq(ancestral_reads, ancestral_fastq_path, 'ANCESTRAL')

    inseq_fastq_path = output_prefix + '.inseq_reads'
    write_reads_to_fastq(inseq_reads, inseq_fastq_path, 'INSEQ')

    combined_fastq_path = output_prefix + '.combined'
    combine_fastq_files(ancestral_fastq_path, inseq_fastq_path, combined_fastq_path)


def get_reads(genome, read_length, fragment_length, step_size):
    reads = []
    positions = [i for i in range(0, len(genome), step_size) if i+fragment_length < len(genome)]

    for i in positions:
        forward_start = i
        forward_end = i+read_length
        reverse_start = i+fragment_length-read_length
        reverse_end = i+fragment_length

        forward_read = genome[forward_start:forward_end]
        reverse_read = revcomp(genome[reverse_start:reverse_end])

        reads.append( (forward_read, reverse_read) )

    return reads


def generate_random_sequence(inseq_length):
    seq = ''
    for i in range(inseq_length):
        seq += int_to_nuc(randint(1,4))
    return seq


def insert_sequence_center(genome, inseq, direct_repeat_length):
    half_length = int(len(genome) / 2)
    inserted = genome[:half_length] + inseq + genome[(half_length-direct_repeat_length):half_length] + genome[half_length:]
    return inserted


def int_to_nuc(integer):
    if integer == 1:
        return 'A'
    elif integer == 2:
        return 'C'
    elif integer == 3:
        return 'G'
    elif integer == 4:
        return 'T'
    else:
        print("WEIRD ERROR")
        sys.exit()

def revcomp(read):
    reversed_seq = ''
    for l in reversed(read):
        if l == 'A':
            reversed_seq += 'T'
        elif l == 'T':
            reversed_seq += 'A'
        elif l == 'C':
            reversed_seq += 'G'
        elif l == 'G':
            reversed_seq += 'C'
        else:
            print("WEIRD ERROR")
            sys.exit()
    return reversed_seq


def write_sequence(sequence, name, out_path):
    record = SeqRecord(Seq(sequence, IUPAC.IUPACUnambiguousDNA), id = name, description=name)

    with open(out_path, "w") as output_handle:
        SeqIO.write([record], output_handle, "fasta")


def write_reads_to_fastq(reads, out_prefix, read_prefix):
    out_forward = open(out_prefix + '.R1.fq', 'w')
    out_reverse = open(out_prefix + '.R2.fq', 'w')
    for i in range(len(reads)):
        fwd_out = '@%s_%d\n' % (read_prefix, i) + reads[i][0]+'\n' + '+\n' + 'K'*len(reads[i][0])+'\n'
        rev_out = '@%s_%d\n' % (read_prefix, i) + reads[i][1] + '\n' + '+\n' + 'K' * len(reads[i][1])+'\n'

        out_forward.write(fwd_out)
        out_reverse.write(rev_out)

    out_forward.close()
    out_reverse.close()

def combine_fastq_files(fastq_prefix1, fastq_prefix2, combined_prefix):
    shell("cat %s.R1.fq %s.R1.fq > %s.R1.fq" % (fastq_prefix1, fastq_prefix2, combined_prefix))
    shell("cat %s.R2.fq %s.R2.fq > %s.R2}.fq" % (fastq_prefix1, fastq_prefix2, combined_prefix))

if __name__ == '__main__':
    main()