#!/usr/bin/python3
# -*- coding: utf-8 -*-
'''
© 2017 Synchon Mandal
'''
import os, argparse, textwrap, csv
from operator import itemgetter
from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord
from collections import OrderedDict
import pandas as pd
import matplotlib.pyplot as plt

# global variables
PATH = os.getcwd()
ACCESSION_ID = ""
CODON_LIST = ()

def perform_search_for_accession_id(handle_id, email_id):
    '''searches the local directory /input for the .fasta file (offline) and if not found falls back to the Entrez database (online)'''

    global PATH

    if os.path.exists(PATH + f'/input/{handle_id}/{handle_id}.fasta'):
        # path change and data extraction (offline)
        os.chdir(PATH + f'/input/{handle_id}')
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

    global CODON_LIST
    # Generates codons and stores it in a list for analysis
    CODON_LIST = compute_codon(dna)

    # Converts DNA to protein
    protein = SeqRecord(dna.translate())

    return protein

def write_protein_sequence_to_file(protein):
    '''writes a fasta file in the /output directory'''

    global PATH
    global ACCESSION_ID

    outputPath = PATH + f'/output/{ACCESSION_ID}'

    # Checks for output directory
    if not os.path.exists(outputPath):
        os.makedirs(outputPath)

    os.chdir(PATH + f'/output/{ACCESSION_ID}')
    SeqIO.write(protein, f'{ACCESSION_ID}_output.fasta','fasta')
    os.chdir(PATH)

def analyse_mutation_dataset(dataset_file):
    '''analyses the given mutation dataset and returns the analysed dataset list'''

    global PATH
    global CODON_LIST

    # input path
    os.chdir(PATH + f'/input/{dataset_file}/')

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
        row.append(CODON_LIST[row[1]])
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
    os.chdir(PATH + f'/output/{dataset_file}/')

    return analysed_dataset

def plot_graph_of_analysed_dataset(dataset_list):
    '''plots the required graphs for the information from the mutation dataset'''

    global ACCESSION_ID
    # change it for different graph style; default: fivethirtyeight
    plt.style.use('fivethirtyeight')

    # only change the default values when needed
    analysed_dataframe = pd.DataFrame(data=dataset_list)
    analysed_dataframe.to_csv(f'{ACCESSION_ID}_output.csv', encoding='utf-8')
    analysed_dataframe.plot(title='Dataset Analysis', grid=True, fontsize=10, figsize=(20,6), x=['Original Nucleotide','Mutated Nucleotide'], y='Occurences', kind='bar', width=0.9, alpha=0.75)
    plt.savefig(f'{ACCESSION_ID}.pdf',bbox_inches='tight')

def main():
    global ACCESSION_ID

    # Command-line argument parser
    parser = argparse.ArgumentParser()
    parser.add_argument('id', type=str, help='NCBI accession number')
    parser.add_argument('email', type=str, help='NCBI registered email-id')
    args = parser.parse_args()

    ACCESSION_ID = args.id

    # Passing the arguments from command line and passing to the required functions
    write_protein_sequence_to_file(convert_DNA_to_protein(perform_search_for_accession_id(handle_id=args.id, email_id=args.email)))
    plot_graph_of_analysed_dataset(analyse_mutation_dataset(dataset_file=args.id))

if __name__ == '__main__':
    main()
