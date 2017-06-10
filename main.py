#!/usr/local/bin/python3.6
'''
© 2017 Synchon Mandal
'''

import sys, os, getopt, xlrd, xlsxwriter, textwrap
from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord
from collections import OrderedDict
import matplotlib.pyplot as plt
import pandas as pd

#global variables
root_path = os.getcwd()
email_id = ""
accession_id = ""
nucleotide_table = ()

def usage():
        print('Usage: ./main.py -a accession-number -e email-id\n')
        print('-h --help               - check the available options\n')
        print('-a --accession_id=id    - NCBI accession number\n')
        print('-e --email_id=email     - NCBI registered email-id\n')
        print('Example:\n')
        print('./main.py -a 123ABC -e abc@efg.syz\n')
        sys.exit(0)

def command_line_argument_check():
    '''checks for command-line arguments and performs necessary opearation'''
    if not len(sys.argv) > 1 :
         usage()

    try:
        opts, args = getopt.getopt(sys.argv[1:],'ha:e:',['help','local_file=','accession_id=','email_id='])
    except getopt.GetoptError as err:
        sys.stderr.write(err)
        usage()

    for o,a in opts:
        if o in ("-h","--help"):
            usage()

        elif o in ("-a","--accession_id"):
            global accession_id
            accession_id = a
            if os.path.exists(root_path + f'/input/{accession_id}'):
                write_protein_sequence_to_file(convert_DNA_to_protein(perform_local_search(file_name=accession_id)))
                mutation_analysis(accession_id=accession_id)
            else:
                write_protein_sequence_to_file(convert_DNA_to_protein(perform_online_Entrez_search(handle_id=accession_id)))
                mutation_analysis(accession_id=accession_id)

        elif o in ("-e","--email_id"):
            email_id = a

        else:
            assert False,"Unhandled Option"

def perform_online_Entrez_search(handle_id):
    '''searches the Entrez database and returns the result'''
    global email_id
    Entrez.email = email_id
    handle = Entrez.efetch(db='nucleotide', rettype='fasta', retmode='text', id =f'{handle_id}')
    seq_record = SeqIO.read(handle, 'fasta')
    handle.close()
    return seq_record.seq

def perform_local_search(file_name):
    '''parses the local file stored in the /input directory'''
    global root_path
    #input path
    os.chdir(root_path + f'/input/{file_name}')
    #data extraction
    seq_record = SeqIO.read(f'{file_name}.fasta','fasta')
    #root
    os.chdir(root_path)
    return seq_record.seq

def calculate_codon(nucleotide):
    '''splits the DNA into codons'''
    return tuple(textwrap.wrap(str(nucleotide),3))

def convert_DNA_to_protein(dna):
    '''translates the DNA to Protein'''
    global nucleotide_table
    nucleotide_table = calculate_codon(dna)
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
    '''analyses the given dataset for mutations and extracts the data for graphical representation'''
    global root_path
    #input path
    os.chdir(root_path + f'/input/{dataset_file}/')
    #opening the dataset and accessing it
    workbook = xlrd.open_workbook(f'{dataset_file}.xlsx')
    worksheet = workbook.sheet_by_index(0)
    #extracting the data
    temporaryList = []
    cellValueList = []
    for row in range(1,worksheet.nrows):
        for col in range(worksheet.ncols):
            temporaryList.append(worksheet.cell(row,col).value)
        cellValueList.append(temporaryList)
        temporaryList = []
    #initiating parameters
    countDictionary = {}
    positionList = []
    positionDictionary = {}
    preCountList = []

    #construction of the primary list for counting
    for value in cellValueList:
        appendedData = value[0] + value[2] + value[3]
        value.append(appendedData)
        preCountList.append(appendedData)

    #construction of set and operation for counting
    preCountSet = set(preCountList)

    for setValue in preCountSet:
        countDictionary[setValue] = preCountList.count(setValue)

    #positon operation
    for setValue in preCountSet:
        for listValue in cellValueList:
            if setValue == listValue[4]:
                positionList.append(int(listValue[1]))
        positionDictionary[setValue] = positionList
        positionList = []

    countDictionary = OrderedDict(sorted(countDictionary.items(), key=lambda t: t[0]))
    positionDictionary = OrderedDict(sorted(positionDictionary.items(), key=lambda t: t[0]))

    #parameters for further analysis
    aminoAcidList = []
    mutatedAminoAcidList = []
    mutatedNucleotideList = []
    occurencesList = []
    positionList = []
    #converts the dictionary items into tuple
    countDictionaryItems = tuple(countDictionary.items())
    positionDictionaryItems = tuple(positionDictionary.items())

    for countItem, positionItem in zip(countDictionaryItems,positionDictionaryItems):
            aminoAcidList.append(countItem[0][0])
            mutatedAminoAcidList.append(countItem[0][1])
            mutatedNucleotideList.append(countItem[0][2:])
            occurencesList.append(countItem[1])
            positionList.append(positionItem[1])

    #changes path to the /output directory; no further change required
    os.chdir(root_path + f'/output/{dataset_file}/')

    return countDictionary, positionDictionary, aminoAcidList, mutatedAminoAcidList, mutatedNucleotideList, occurencesList, positionList


def render_worksheet_of_analysed_dataset(count_dictionary, position_dictionary):
    global nucleotide_table
    global accession_id

    #worksheet environment setup
    workbook = xlsxwriter.Workbook(f'{accession_id}.xlsx')
    worksheet = workbook.add_worksheet('Mutation Dataset Analysis')

    #defines parameters
    heading = workbook.add_format({'font_name': 'Courier New', 'bold': True, 'align': 'center', 'valign': 'vcenter', 'text_wrap': True})
    values = workbook.add_format({'font_name': 'Courier New', 'align': 'center', 'valign': 'vcenter'})
    worksheet.set_column(0, 5, 30)

    #worksheet operation
    worksheet.write(0, 0, 'Amino Acid', heading)
    worksheet.write(0, 1, 'Mutated Amino Acid', heading)
    worksheet.write(0, 2, 'Nucleotide', heading)
    worksheet.write(0, 3, 'Mutated Nucleotide', heading)
    worksheet.write(0, 4, 'Occurences', heading)
    worksheet.write(0, 5, 'Positions', heading)

    #fills the data
    rowForFillingCount = 0
    for key, value in count_dictionary.items():
        rowForFillingCount += 1
        worksheet.write(rowForFillingCount, 0, key[0], values)
        worksheet.write(rowForFillingCount, 1, key[1], values)
        worksheet.write(rowForFillingCount, 3, key[2:], values)
        worksheet.write(rowForFillingCount, 4, value, values)


    rowForFillingPosition = 0
    nucleotideList = []
    tempNucleotideList = []
    for key, value in position_dictionary.items():
        rowForFillingPosition += 1
        for individualPosition in value:
            tempNucleotideList.append(nucleotide_table[int(individualPosition)])
            nucleotideString = ','.join(str(codon) for codon in tempNucleotideList)

        nucleotideList.append(nucleotideString)
        worksheet.write(rowForFillingPosition, 2, nucleotideString, values)
        tempNucleotideList = []

        value = ','.join(str(number) for number in value)
        worksheet.write(rowForFillingPosition, 5, value, values)

    #closing the file
    workbook.close()

    return nucleotideList

def plot_graphs_of_analysed_dataset(amino_acid_list, mutated_amino_acid_list, nucleotide_list, mutated_nucleotide_list, occurences_list, position_list):
    '''plots the required graphs for the information from the mutation dataset'''
    global accession_id
    #change it for different graph style; default: fivethirtyeight
    plt.style.use('fivethirtyeight')

    #AMINO ACID; only change the default values when needed
    aminoAcidAnalysisDataFrame = pd.DataFrame({'Amino Acid': amino_acid_list, 'Mutated Amino Acid': mutated_amino_acid_list, 'Occurences': occurences_list, 'Positions': position_list})
    aminoAcidAnalysisDataFrame.plot(title='Amino Acid Dataset Analysis', grid=True, fontsize=10, figsize=(20,6), x=['Amino Acid','Mutated Amino Acid'], y='Occurences', kind='bar', width=0.9, alpha=0.75)
    plt.savefig(f'{accession_id}_amino_acid.pdf',bbox_inches='tight')

    #NUCLEOTIDE; only change the default values when needed
    nucleotideAnalysisDataFrame = pd.DataFrame({'Nucleotide': nucleotide_list, 'Mutated Nucleotide': mutated_nucleotide_list, 'Occurences': occurences_list, 'Positions': position_list})
    nucleotideAnalysisDataFrame.plot(title='Nucleotide Dataset Analysis', grid=True, fontsize=10, figsize=(20,6), x=['Nucleotide','Mutated Nucleotide'], y='Occurences', kind='bar', width=0.9, alpha=0.75)
    plt.savefig(f'{accession_id}_nucleotide.pdf',bbox_inches='tight')

def mutation_analysis(accession_id):
    '''it analyses the mutation dataset provided using the helper functions'''
    counting_dictionary, position_dictionary, amino_acid, mutated_amino_acid, mutated_nucleotide, occurences, position = analyse_mutation_dataset(dataset_file=accession_id)
    nucleotide = render_worksheet_of_analysed_dataset(count_dictionary=counting_dictionary,position_dictionary=position_dictionary)
    plot_graphs_of_analysed_dataset(amino_acid_list=amino_acid,mutated_amino_acid_list=mutated_amino_acid,nucleotide_list=nucleotide,mutated_nucleotide_list=mutated_nucleotide,occurences_list=occurences,position_list=position)

def main():
    command_line_argument_check()

if __name__ == '__main__':
    main()
