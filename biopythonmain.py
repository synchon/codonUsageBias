"""
Synchon Mandal
2017
"""
import sys
from operator import lt
from Bio import Entrez, SeqIO

def usage():
        sys.stderr.write("Usage: python3 main.py <email-id> <accession-no>\n")
        sys.exit(0)

def main():

    #Entrez search
    def onlineEntrezSearch(for= str(sys.argv[2]),by= str(sys.argv[1])):
        Entrez.email = str(sys.argv[1])
        handle_id = str(sys.argv[2])
        handle = Entrez.efetch(db="nucleotide", rettype="fasta", retmode="text", id=handle_id)
        seq_record = SeqIO.read(handle, "fasta")
        handle.close()

    #local parsing
    def Local
    #DNA and Protein operations
    dna = seq_record.seq
    protein = dna.translate()
    print(f'No. of codons: {len(protein)}')
    print(dna)
    print(protein)


if __name__ == '__main__':
    #command-line argument check functions
    def argumentCheckFailed():
        usage()

    def argumentCheckPassed():
        main()

    (lt(len(sys.argv),3) and argumentCheckFailed()) or (argumentCheckPassed())
