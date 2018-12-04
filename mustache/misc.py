import warnings
warnings.filterwarnings("ignore")
import sys
from scipy.stats import poisson


def get_bam_contig_dict(bam_file):
    contig_dict = {}
    for contig in bam_file.header['SQ']:
        contig_dict[contig['SN']] = contig['LN']
    return contig_dict

def revcomp(read):
    reversed_seq = ''
    for l in reversed(read.upper()):
        if l == 'A':
            reversed_seq += 'T'
        elif l == 'T':
            reversed_seq += 'A'
        elif l == 'C':
            reversed_seq += 'G'
        elif l == 'G':
            reversed_seq += 'C'
        else:
            reversed_seq += l
    return reversed_seq

def poisson_test_greater(x, mu):

    return 1 - poisson.cdf(x-1, mu)


if __name__ == "__main__":
    print(poisson_test_greater(1, mu=1))