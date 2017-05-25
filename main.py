#!/usr/local/bin/python3.6
'''
Synchon Mandal
2017
'''

import sys, os, getopt
from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord

#global variables
RootPath = os.getcwd()
email_id = ""
accession_id = ""

def usage():
        print('Usage: ./main.py -a accession-number -e email-id\n')
        print('-h --help               - check the available options\n')
        print('-a --accession_id=id    - NCBI accession number\n')
        print('-e --email_id=email     - NCBI registered email-id\n')
        print('Example:\n')
        print('python3 main.py -a 123ABC -e abc@efg.syz\n')
        sys.exit(0)

def commandLineArgumentCheck():
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
            if os.path.exists(RootPath + f'input/{accession_id}'):
                writeProteinSequenceToFile(convertDNAtoProtein(localFileParsing(file_name=accession_id)))
            else:
                writeProteinSequenceToFile(convertDNAtoProtein(onlineEntrezSearch(handle_id=accession_id)))
        elif o in ("-e","--email_id"):
            email_id = a
        else:
            assert False,"Unhandled Option"

def onlineEntrezSearch(handle_id):
    '''searches the Entrez database and returns the result'''
    global email_id
    Entrez.email = email_id
    handle = Entrez.efetch(db='nucleotide', rettype='fasta', retmode='text', id =f'{handle_id}')
    seq_record = SeqIO.read(handle, 'fasta')
    handle.close()
    return seq_record.seq

def localFileParsing(file_name):
    '''parses the local file stored in the /input directory'''
    global RootPath
    os.chdir(RootPath + f'/input/{file_name}')
    seq_record = SeqIO.read(f'{file_name}.fasta','fasta')
    os.chdir(RootPath)
    return seq_record.seq

def convertDNAtoProtein(dna):
    '''translates the DNA to Protein'''
    protein = SeqRecord(dna.translate())
    # print(f'No. of codons: {len(protein)}')
    return protein

def writeProteinSequenceToFile(protein):
    '''writes a fasta file in the /output directory'''
    global RootPath
    global accession_id
    outputPath = RootPath + f'/output/{accession_id}'
    if not os.path.exists(outputPath):
        os.makedirs(outputPath)
    os.chdir(RootPath + f'/output/{accession_id}')
    seq_record = SeqIO.write(protein,f'{accession_id.lower()}_output.fasta','fasta')
    os.chdir(RootPath)

def main():
    commandLineArgumentCheck()

if __name__ == '__main__':
    main()
