#!/usr/bin/python3
# -*- coding: utf-8 -*-
'''
© 2017 Synchon Mandal
'''

import sys, os, getopt, textwrap, csv
from operator import itemgetter
from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord
from collections import OrderedDict
import pandas as pd
import matplotlib.pyplot as plt

# global variables
root_path = os.getcwd()
email_id = ""
accession_id = ""
codon_list = ()

def usage():
        print('Usage: ./main.py -a accession-number -e email-id\n')
        print('-h --help               - check the available options\n')
        print('-a --accession_id=id    - NCBI accession number\n')
        print('-e --email_id=email     - NCBI registered email-id\n')
        print('Example:\n')
        print('./main.py -a 123ABC -e user@example.com\n')
        sys.exit(0)

def command_line_argument_check():
    '''checks for command-line arguments and performs necessary operation'''
    if not len(sys.argv) > 1 :
         usage()

    try:
        opts, args = getopt.getopt(sys.argv[1:],'ha:e:',['help','accession_id=','email_id='])
    except getopt.GetoptError as err:
        sys.stderr.write(err)
        usage()

    for o,a in opts:
        if o in ("-h","--help"):
            usage()

        elif o in ("-a","--accession_id"):
            global accession_id
            accession_id = a
            write_protein_sequence_to_file(convert_DNA_to_protein(perform_search_for_accession_id(handle_id=accession_id)))
            plot_graph_of_analysed_dataset(analyse_mutation_dataset(dataset_file=accession_id))

        elif o in ("-e","--email_id"):
            email_id = a

        else:
            assert False,"Unhandled Option"

def perform_search_for_accession_id(handle_id):
    '''searches the local directory /input for the .fasta file (offline) and if not found falls back to the Entrez database (online)'''
    global root_path
    global email_id

    if os.path.exists(root_path + f'/input/{accession_id}/{accession_id}.fasta'):
        # path change and data extraction (offline)
        os.chdir(root_path + f'/input/{handle_id}')
        seq_record = SeqIO.read(f'{handle_id}.fasta','fasta')
    else:
        # based on Entrez API (online)
        Entrez.email = email_id
        handle = Entrez.efetch(db='nucleotide', rettype='fasta', retmode='text', id =f'{handle_id}')
        seq_record = SeqIO.read(handle, 'fasta')
        handle.close()

    return seq_record.seq

def compute_codon(nucleotide):
    '''splits the DNA sequence into codons'''
    return tuple(textwrap.wrap(str(nucleotide),3))
    # ['AGAGAG...'] => ['AGA','GAG','...',...]

def convert_DNA_to_protein(dna):
    '''translates the DNA to Protein'''
    global codon_list
    codon_list = compute_codon(dna)
    protein = SeqRecord(dna.translate())
    # print(f'No. of codons: {len(protein)}')
    return protein

def write_protein_sequence_to_file(protein):
    '''writes a fasta file in the /output directory'''
    global root_path
    global accession_id
    outputPath = root_path + f'/output/{accession_id}'
    if not os.path.exists(outputPath):
        os.makedirs(outputPath)
    os.chdir(root_path + f'/output/{accession_id}')
    seq_record = SeqIO.write(protein,f'{accession_id.lower()}_output.fasta','fasta')
    os.chdir(root_path)

def analyse_mutation_dataset(dataset_file):
    '''analyses the given mutation dataset and returns the analysed dataset list'''
    global root_path
    global codon_list
    # input path
    os.chdir(root_path + f'/input/{dataset_file}/')
    # CSV parsing
    with open(f'{dataset_file}.csv','r') as mutation_file:
        mutation_reader = csv.reader(mutation_file)
        next(mutation_reader)
        mutation_dataset = list(mutation_reader)
        # [['K', '3', 'E', 'GAG'],...]

    # initiating parameters
    original_and_mutated_nucleotide_list = []
    analysed_dataset = []

    # row-based operation on dataset
    for row in mutation_dataset:
        row[1] = int(row[1])
        row.append(codon_list[row[1]])
        original_and_mutated_nucleotide = row[4] + row[3]
        row.append(original_and_mutated_nucleotide)
        original_and_mutated_nucleotide_list.append(original_and_mutated_nucleotide)
    # ['A','3','G','GAG'] => ['A',3,'G','GAG','AAG','AAGGAG']

    # construction of set for removing duplicates
    original_and_mutated_nucleotide_set = set(original_and_mutated_nucleotide_list)
    # ['AAGGAG','AAGGAG',...] => {'AAGGAG',...}

    for set_value in original_and_mutated_nucleotide_set:
        position_list = []
        occurences = original_and_mutated_nucleotide_list.count(set_value)
        original_nucleotide = set_value[:3]
        mutated_nucleotide = set_value[3:]

        for row in mutation_dataset:
            if set_value == row[5]:
                position_list.append(row[1])
                original_amino_acid = row[0]
                mutated_amino_acid = row[2]

        final_data = OrderedDict()
        final_data['Original Nucleotide'] = original_nucleotide
        final_data['Mutated Nucleotide'] = mutated_nucleotide
        final_data['Occurences'] = occurences
        final_data['Position'] = position_list
        final_data['Original Amino Acid'] = original_amino_acid
        final_data['Mutated Amino Acid'] = mutated_amino_acid
        analysed_dataset.append(final_data)
        analysed_dataset = sorted(analysed_dataset, key=itemgetter('Original Nucleotide'))
        # [OrderedDict([('Original Nucleotide', 'AAC'),
        #               ('Mutated Nucleotide', 'GAC'),
        #               ('Occurences', 1),
        #               ('Position', [139]),
        #               ('Original Amino Acid', 'N'),
        #               ('Mutated Amino Acid', 'D')]),
        #  OrderedDict([('Original Nucleotide', 'AAC'),
        #               ('Mutated Nucleotide', 'CAC'),
        #               ('Occurences', 1),
        #               ('Position', [139]),
        #               ('Original Amino Acid', 'N'),
        #               ('Mutated Amino Acid', 'H')]),...]

    # changes path to the /output directory; no further directory change required
    os.chdir(root_path + f'/output/{dataset_file}/')

    return analysed_dataset

def plot_graph_of_analysed_dataset(dataset_list):
    '''plots the required graphs for the information from the mutation dataset'''
    global accession_id
    # change it for different graph style; default: fivethirtyeight
    plt.style.use('fivethirtyeight')

    # only change the default values when needed
    analysed_dataframe = pd.DataFrame(data=dataset_list)
    analysed_dataframe.to_csv(f'{accession_id}_output.csv', encoding='utf-8')
    analysed_dataframe.plot(title='Dataset Analysis', grid=True, fontsize=10, figsize=(20,6), x=['Original Nucleotide','Mutated Nucleotide'], y='Occurences', kind='bar', width=0.9, alpha=0.75)
    plt.savefig(f'{accession_id}.pdf',bbox_inches='tight')

def main():
    command_line_argument_check()

if __name__ == '__main__':
    main()
