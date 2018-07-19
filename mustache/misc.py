import sys


def get_bam_contig_dict(bam_file):
    contig_dict = {}
    for contig in bam_file.header['SQ']:
        contig_dict[contig['SN']] = contig['LN']
    return contig_dict

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
            reversed_seq += l
    return reversed_seq
