"""
2017
Synchon Mandal
"""
#!/usr/local/bin/python3.6

import os, sys, datetime, time, computation
#from multiprocessing import Pool

def usage():
        sys.stderr.write("Usage: python3 main.py <input_file> <mutation: mut/no-mut>\n")
        sys.exit(0)

def main():

    #defines the root directory path
    rootPath = os.getcwd()
    #sets the output file name
    OutputFileName = InputFileName + str(datetime.datetime.now())
    os.makedirs(rootPath + f'/output/{OutputFileName}')

    #CODON GENERATION
    #runs the sorting function for codon conversion
    print('\nSorting the genome file...')
    StartTimeOfSorting = time.time()
    sortedTuple = computation.sortTextFile(fileName=IutputFileName)
    print (f'\nTotal number of codons: {len(sortedTuple)}')
    totalTimeOfSorting: float = round(time.time() - StartTimeOfSorting, 3)
    print(f'\nSorted.\nTotal time required for sorting: {totalTimeOfSorting} seconds.\n==========================================')


    #CONVERSION TO AMINO ACID
    #runs the computation function for the data obtained
    print('\nTranslating into amino acids...')
    startTimeOfTranslating: float = time.time()
    oneLetterTuple: Tuple = computation.makeCodonFromSortedFile(codonTuple=sortedTuple)
    totalTimeOfTranslating: float = round(time.time() - startTimeOfTranslating, 3)
    print(f'\nTranslated.\nTotal time required for translating: {totalTimeOfTranslating} seconds.\n==========================================')


    #TXT GENERATION
    #writes the data into a .txt file
    computation.renderTextFile(fileName=OutputFileName,aminoAcidSequence=oneLetterTuple)
    print('\n.txt file generated! Check the output/txt directory.\n==========================================\nMutation Analysis')


    #MUTATION
    #muatation computation
    if responseForMutation == 'no-mut':
        pass

    elif responseForMutation == 'mut':
        #prints out the data obtained for further processing
        #TABLE
        computation.viewTableOfAminoAcid(aminoAcidSequence=oneLetterTuple)

        #accepts the position number(s) and amino acid(s)
        positionNumbers: Tuple = tuple(input('\nPlease provide the position number(s) separated by commas:\n').split(','))
        positionRequiredAminoAcids: Tuple = tuple(input('\nPlease enter the required amino acid(s) separated by commas:\n').split(','))
        print('\n==========================================')
        #returns the mutation results
        positionCodon, positionAminoAcid, finalMutatedAminoAcidDictionary = computation.userDefinedMutation(DNASequence=sortedTuple,aminoAcidSequence=oneLetterTuple,positionNumbersOfAminoAcids=positionNumbers,positionRequiredMutatedAminoAcids=positionRequiredAminoAcids)

        #error checking
        if positionCodon == None or positionAminoAcid == None or finalMutatedAminoAcidDictionary == None:
            return
        else:
            #WORKSHEET GENERATION
            #renders the worksheet file
            computation.renderWorksheetOfUserDefinedMutation(fileName=OutputFileName,DNAdictionary=positionCodon,AminoAcidDictonary=positionAminoAcid,MutatedAminoAcidDictionary=finalMutatedAminoAcidDictionary)
            print('\nWorksheet generated! Check the output/userDefinedMutationWorksheet directory.\n==========================================')

    #DATASET ANALYSIS
    #analyses, operates and plots the Dataset file
    dictOfCount, dictOfPosition, aminoAcid, mutatedAminoAcid, mutatedNucleotide, occurences, position = computation.analysisForMutationDataset(fileName=IutputFileName)
    nucleotide = computation.renderWorksheetForMutationDatasetAnalysis(fileName=IutputFileName,outputxlsxFile=OutputFileName,countDictionary=dictOfCount,positionDictionary=dictOfPosition)
    print('\nDataset Analysis worksheet generated! Check the output/datasetWorksheet directory.\n==========================================')
    computation.plotGraphForMutationDatasetAnalysis(fileName=OutputFileName,aminoAcidList=aminoAcid,mutatedAminoAcidList=mutatedAminoAcid,nucleotideAnalysisList=nucleotide,mutatedNucleotideAnalysisList=mutatedNucleotide,occurencesList=occurences,positionList=position)
    print(f'\nDataset Analysis graph generated! Check the output/datasetAnalysisGraph/{OutputFileName} directory.\nDone.\n==========================================')

if __name__ == '__main__':
    #command-line argument check
    if len(sys.argv) < 3:
        usage()
    else:
        IutputFileName: str = str(sys.argv[1]).upper()
        responseForMutation: str = str(sys.argv[2]).lower()
        main()
