import sys
from snakemake import shell
import pandas as pd

def blast_fasta(fasta, blastdb, outfile, threads=2):
    add_header = 'echo -e "qseqid\tsseqid\tstitle\tpident\tlength\tmismatch\tgapopen' \
                 '\tqstart\tqend\tsstart\tsend\tevalue\tbitscore" > {outfile}'.format(outfile=outfile)
    shell(add_header)

    blast_command = 'blastn -query {fasta} -db {blastdb} -evalue 0.01 -outfmt "6 ' \
                    'qseqid sseqid stitle pident length mismatch gapopen qstart qend ' \
                    'sstart send evalue bitscore" -num_threads {threads} >> {outfile}'.format(
        fasta=fasta, threads=threads, outfile=outfile, blastdb=blastdb
    )
    shell(blast_command)

def process_blast_results(blast_results, evalue_cutoff=0.01):

    results = dict()
    with open(blast_results) as infile:
        header = infile.readline().strip().split('\t')
        for l in infile:
            line = l.strip().split('\t')
            line = {header[i]:line[i] for i in range(len(line))}
            qseqid = line['qseqid']
            stitle = line['stitle']
            evalue = float(line['evalue'])
            is_name = stitle.split()[1]
            is_family = stitle.split()[2]
            index = int(qseqid.split('_')[0])
            orient = qseqid.split('_')[-1]

            if index in results and orient in results[index]:
                current_evalue = results[index][orient]['evalue']
                if evalue > current_evalue:
                    results[index][orient]['evalue'] = evalue
                    results[index][orient]['is_name'] = set([is_name])
                    results[index][orient]['is_family'] = set([is_family])
                elif evalue == current_evalue:
                    results[index][orient]['is_name'].add(is_name)
                    results[index][orient]['is_family'].add(is_family)
            elif index in results:
                results[index][orient] = dict()
                results[index][orient]['evalue'] = evalue
                results[index][orient]['is_name'] = set([is_name])
                results[index][orient]['is_family'] = set([is_family])
            else:
                results[index] = dict()
                results[index][orient] = dict()
                results[index][orient]['evalue'] = evalue
                results[index][orient]['is_name'] = set([is_name])
                results[index][orient]['is_family'] = set([is_family])

    for index in results:
        for orient in results[index]:
            results[index][orient]['is_name'] = ';'.join(list(results[index][orient]['is_name']))
            results[index][orient]['is_family'] = ';'.join(list(results[index][orient]['is_family']))

    header = ['evalue_5p', 'evalue_3p', 'IS_5p', 'IS_3p', 'IS_family_5p', 'IS_family_3p']
    final_results = dict()
    for index in results:
        if '5p' in results[index] and '3p' in results[index]:
            final_results[index] = [results[index]['5p']['evalue'], results[index]['3p']['evalue'],
                                    results[index]['5p']['is_name'], results[index]['3p']['is_name'],
                                    results[index]['5p']['is_family'], results[index]['3p']['is_family']]
        elif '5p' in results[index]:
            final_results[index] = [results[index]['5p']['evalue'], None,
                                    results[index]['5p']['is_name'], None,
                                    results[index]['5p']['is_family'], None]
        elif '3p' in results[index]:
            final_results[index] = [None, results[index]['3p']['evalue'],
                                    None, results[index]['3p']['is_name'],
                                    None, results[index]['3p']['is_family']]

    final_results = pd.DataFrame.from_dict(final_results , orient='index', columns=header)

    return final_results
