import warnings
warnings.filterwarnings("ignore")
import sys
from snakemake import shell
import pandas as pd
from random import randint
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from os.path import isfile
from os.path import join


# conda install -c biocore hmmer

def run_hmmsearch(fasta, outfile, database, threads=2):
    hmmsearch_command = 'hmmsearch --tblout={outfile} --max --noali --cpu={threads} {database} {fasta} 1> /dev/null'.format(
        fasta=fasta, outfile=outfile, threads=threads, database=database
    )
    shell(hmmsearch_command)

# Figure out evalue_cutoff
def process_hmm_results(hmm_results, evalue_cutoff=0.01):

    results = dict()

    if isfile(hmm_results):

        with open(hmm_results) as infile:
            for i in range(3):
                skip = infile.readline()
            header = ["target", "accession", "query", "accession", "E-value", "score", "bias"]
            for l in infile:
                if "#" not in l:
                    line = l.strip().split()[0:7]
                    line = {header[i]: line[i] for i in range(len(line))}

                    identifier = line['target']
                    is_name = line['query'].split(".")[0]
                    evalue = float(line['E-value'])

                    if evalue > evalue_cutoff:
                        continue

                    # strand = identifier.split("_")[-1] # + or -
                    orient = identifier.split("_")[1]  # 5p or 3p
                    index = int(identifier.split("_")[0])
                    start_pos = identifier.split("_")[2]
                    end_pos = identifier.split("_")[3]
                    if index in results and orient in results[index]:
                        current_evalue = results[index][orient]['evalue']
                        if evalue < current_evalue:
                            results[index][orient]['evalue'] = evalue
                            results[index][orient]['is_name'] = set([is_name])
                            results[index][orient]['start_pos'] = set([start_pos])
                            results[index][orient]['end_pos'] = set([end_pos])
                        elif evalue == current_evalue:
                            results[index][orient]['is_name'].add(is_name)
                    elif index in results:
                        results[index][orient] = dict()
                        results[index][orient]['evalue'] = evalue
                        results[index][orient]['is_name'] = set([is_name])
                        results[index][orient]['start_pos'] = set([start_pos])
                        results[index][orient]['end_pos'] = set([end_pos])
                    else:
                        results[index] = dict()
                        results[index][orient] = dict()
                        results[index][orient]['evalue'] = evalue
                        results[index][orient]['is_name'] = set([is_name])
                        results[index][orient]['start_pos'] = set([start_pos])
                        results[index][orient]['end_pos'] = set([end_pos])

    for index in results:
        for orient in results[index]:
            results[index][orient]['is_name'] = ';'.join(list(results[index][orient]['is_name']))
            results[index][orient]['start_pos'] = ';'.join(list(results[index][orient]['start_pos']))
            results[index][orient]['end_pos'] = ';'.join(list(results[index][orient]['end_pos']))

        # print(index, orient, results[index][orient]['is_name'], results[index][orient]['evalue'])

    header = ['hmmsearch_evalue_5p', 'hmmsearch_evalue_3p', 'hmmsearch_IS_5p', 'hmmsearch_IS_3p',
              'hmmsearch_start_pos_5p', 'hmmsearch_end_pos_5p', 'hmmsearch_start_pos_3p', 'hmmsearch_end_pos_3p']

    final_results = dict()
    for index in results:
        if '5p' in results[index] and '3p' in results[index]:
            final_results[index] = [results[index]['5p']['evalue'], results[index]['3p']['evalue'],
                                    results[index]['5p']['is_name'], results[index]['3p']['is_name'],
                                    results[index]['5p']['start_pos'], results[index]['5p']['end_pos'],
                                    results[index]['3p']['start_pos'], results[index]['3p']['end_pos']]
        elif '5p' in results[index]:
            final_results[index] = [results[index]['5p']['evalue'], None,
                                    results[index]['5p']['is_name'], None,
                                    results[index]['5p']['start_pos'], results[index]['5p']['end_pos'],
                                    None, None]
        elif '3p' in results[index]:
            final_results[index] = [None, results[index]['3p']['evalue'],
                                    None, results[index]['3p']['is_name'],
                                    None, None,
                                    results[index]['3p']['start_pos'], results[index]['3p']['end_pos']]

    final_results = pd.DataFrame.from_dict(final_results, orient='index', columns=header)
    return final_results
