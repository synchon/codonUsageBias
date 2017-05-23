"""
2017
Synchon Mandal
"""
from typing import List, Dict, Tuple
import json, math, os, xlsxwriter, xlrd
from tabulate import tabulate
from collections import OrderedDict
import matplotlib.pyplot as plt
import pandas as pd

RootPath = os.getcwd()

def CodonTable(fileName):
   #opens the json file for importing Codon Table
   with open(f'{fileName}.json' ,'r') as table:
      data = json.load(table)
      return data

#this is currently for FASTA format only
def ConvertCodon(fileName):

   global RootPath
   #changes to the required path
   os.chdir(RootPath + f'/input/{fileName}/')

   with open(f'{fileName}.txt','rU') as NucleotideFile:
      #skips the first line of FASTA
      next(NucleotideFile)
      #defines the required parameter
      sortedNucleotideList = []
      iterVar = 0
      #sets up the file
      ReadNucleotide = NucleotideFile.read().split('\n')
      ReadNucleotide = ''.join(readNucleotide)

      #ALGORITHM FOR SORTING
      if len(readNucleotide)%3 == 0:
         for codon in range(int(len(readNucleotide)/3)):
            triplet = readNucleotide[iterVar] + readNucleotide[iterVar + 1] + readNucleotide[iterVar + 2]
            sortedNucleotideList.append(triplet)
            iterVar += 3

      elif len(readNucleotide)%3 != 0:
         for codon in range(int(math.floor(len(readNucleotide)/3))):
            try:
                triplet = readNucleotide[iterVar] + readNucleotide[iterVar + 1] + readNucleotide[iterVar + 2]
                sortedNucleotideList.append(triplet)
                iterVar += 3
            except:
                pass
   #changes to the root directory path
   os.chdir(RootPath)
   #converts the list to a tuple
   sortedNucleotideList = tuple(sortedNucleotideList)
   return sortedNucleotideList

def makeCodonFromSortedFile(codonTuple: Tuple) -> Tuple:
   #imports the Codon Table
   codonTable: Dict[str,str] = openJson(fileName='dnaCodonTable')
   #initiates the list
   oneLetterSequence: List[str] = []
   #runs the loop for adding data
   for codon in codonTuple:
      for knownCodonKey, knownCodonValue in codonTable.items():
         if codon == knownCodonKey:
            oneLetterSequence.append(knownCodonValue)
   #converts the list to tuple
   oneLetterSequence = tuple(oneLetterSequence)
   return oneLetterSequence


def renderTextFile(fileName: str, aminoAcidSequence: Tuple) -> None:

    global RootPath
    #sets the path
    os.chdir(RootPath + f'/output/{fileName}/')
    #converts the tuple into a string
    aminoAcidSequence = ''.join(aminoAcidSequence)
    #writes to the file
    with open(f'{fileName}.txt','w',encoding = 'utf-8') as outputFile:
        outputFile.write(aminoAcidSequence)
    #changes to the root directory path
    os.chdir(RootPath)


def viewTableOfAminoAcid(aminoAcidSequence: Tuple) -> None:
    #sets the required parameters
    oneLetterAminoAcidIndex: List[int] = []
    temporaryList: [str,int] = []
    aminoAcidTable: List[[str,int]] = []
    headers: List[str] = ["Position","Amino Acid"]
    #stores the indices in a list
    for indexValue in range(len(aminoAcidSequence)):
        oneLetterAminoAcidIndex.append(indexValue)
    #makes the list required for the table
    for index, value in zip(oneLetterAminoAcidIndex, aminoAcidSequence):
        temporaryList.append(index)
        temporaryList.append(value)
        aminoAcidTable.append(temporaryList)
        temporaryList = []

    print(tabulate(aminoAcidTable, headers, tablefmt="psql"))


def userDefinedMutation(DNASequence: Tuple, aminoAcidSequence: Tuple, positionNumbersOfAminoAcids: Tuple, positionRequiredMutatedAminoAcids: Tuple) -> (Dict[int,int],Dict[int,int],Dict[str,str]):

    #imports codon table
    codonTab: Dict[str,str] = openJson('dnaCodonTable')

    #primary error checking
    if len(positionNumbersOfAminoAcids) == len(positionRequiredMutatedAminoAcids):
        #secondary error checking
        for num in positionNumbersOfAminoAcids:
            if int(num) > len(DNASequence):
                print(f'The index : {num} is not valid!')
                return None, None, None

            else:
                #sets the required parameters
                positionCodon = {}
                positionAminoAcid = {}
                finalMutatedAminoAcidDictionary = {}
                tempAminoAcidList = []
                tempAminoAcidString = ''

                try:
                    #DNA and AMINO ACID
                    for num in positionNumbersOfAminoAcids:
                        positionCodon[int(num)] = DNASequence[int(num)]
                        positionAminoAcid[int(num)] = aminoAcidSequence[int(num)]
                    #MUTATED AMINO ACID
                    for aminoAcid in positionRequiredMutatedAminoAcids:
                        for key, item in codonTab.items():
                            if aminoAcid == item:
                                tempAminoAcidList.append(key)

                        tempAminoAcidString = ','.join(tempAminoAcidList)
                        finalMutatedAminoAcidDictionary[aminoAcid] = tempAminoAcidString
                        tempAminoAcidList = []
                        tempAminoAcidString = ''

                except:
                    pass

                return positionCodon, positionAminoAcid, finalMutatedAminoAcidDictionary

    else:
        print('The number of mutation(s) don\'t match!')
        return None, None, None



def renderWorksheetOfUserDefinedMutation(fileName: str, DNAdictionary: Dict[int,int], AminoAcidDictonary: Dict[int,int], MutatedAminoAcidDictionary: Dict[str,str]) -> None:

    global RootPath
    #sets the required path
    os.chdir(RootPath + f'/output/{fileName}/')

    #worksheet environment setup
    workbook = xlsxwriter.Workbook(f'{fileName}userDefinedMutationWorksheet.xlsx')
    worksheet = workbook.add_worksheet('Mutation Sheet')


    #defines parameters
    heading = workbook.add_format({'font_name': 'Courier New', 'bold': True, 'align': 'center', 'valign': 'vcenter', 'text_wrap': True})
    values = workbook.add_format({'font_name': 'Courier New', 'align': 'center', 'valign': 'vcenter'})
    worksheet.set_column(0, 0, 10)
    worksheet.set_column(1, 3, 20)
    worksheet.set_column(4, 4, 30)
    #worksheet operation
    worksheet.write(0, 0, 'Position', heading)
    worksheet.write(0, 1, 'Nucleotide', heading)
    worksheet.write(0, 2, 'Amino Acid', heading)
    worksheet.write(0, 3, 'Required Amino Acid', heading)
    worksheet.write(0, 4, 'Possible Nucleotides', heading)

    #fills the data
    row = 0
    for key in DNAdictionary:
        row += 1
        worksheet.write(row, 0, key, values)
        worksheet.write(row, 1, DNAdictionary[key], values)

    row = 0
    for key in AminoAcidDictonary:
        row += 1
        worksheet.write(row, 2, AminoAcidDictonary[key], values)

    row = 0
    for key in MutatedAminoAcidDictionary:
        row += 1
        worksheet.write(row, 3, key, values)
        worksheet.write(row, 4, MutatedAminoAcidDictionary[key], values)

    #closes the file
    workbook.close()
    os.chdir(RootPath)


def analysisForMutationDataset(fileName: str) -> (Dict[str,int],Dict[str, List[int]],List[str],List[str],List[str],List[int],List[int]):

    global RootPath
    #sets the path
    os.chdir(RootPath + f'/input/{fileName}/')
    #opening the file and accessing it
    workbook = xlrd.open_workbook(f'{fileName}.xlsx')
    worksheet = workbook.sheet_by_index(0)
    #extracting the data
    temporaryList: List[str] = []
    cellValueList: List[str] = []
    for row in range(1,worksheet.nrows):
        for col in range(worksheet.ncols):
            temporaryList.append(worksheet.cell(row,col).value)
        cellValueList.append(temporaryList)
        temporaryList = []
    #initiating parameters
    countDictionary: Dict[str,int] = {}
    positionList: List[int] = []
    positionDictionary: Dict[str, List[int]] = {}
    preCountList: List[str] = []

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

    os.chdir(RootPath)

    countDictionary = OrderedDict(sorted(countDictionary.items(), key=lambda t: t[0]))
    positionDictionary = OrderedDict(sorted(positionDictionary.items(), key=lambda t: t[0]))

    #parameters for further analysis
    aminoAcidAnalysisList: List[str] = []
    mutatedAminoAcidAnalysisList: List[str] = []
    mutatedNucleotideAnalysisList: List[str] = []
    occurencesAnalysisList: List[int] = []
    positionAnalysisList: List[List[int]] = []
    #converts the dictionary items into tuple
    countDictionaryItems = tuple(countDictionary.items())
    positionDictionaryItems = tuple(positionDictionary.items())

    for countItem, positionItem in zip(countDictionaryItems,positionDictionaryItems):
            aminoAcidAnalysisList.append(countItem[0][0])
            mutatedAminoAcidAnalysisList.append(countItem[0][1])
            mutatedNucleotideAnalysisList.append(countItem[0][2:])
            occurencesAnalysisList.append(countItem[1])
            positionAnalysisList.append(positionItem[1])

    return countDictionary, positionDictionary, aminoAcidAnalysisList, mutatedAminoAcidAnalysisList, mutatedNucleotideAnalysisList, occurencesAnalysisList, positionAnalysisList


def renderWorksheetForMutationDatasetAnalysis(fileName: str, outputxlsxFile, countDictionary, positionDictionary) -> List[str]:
    #imports the required file for operation
    dnaTable: Tuple = sortTextFile(fileName=fileName)

    global RootPath
    #sets the required path
    os.chdir(RootPath + f'/output/{outputxlsxFile}/')
    #worksheet environment setup
    workbook = xlsxwriter.Workbook(f'{outputxlsxFile}.xlsx')
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
    row: int = 0
    for key, value in countDictionary.items():
        row += 1
        worksheet.write(row, 0, key[0], values)
        worksheet.write(row, 1, key[1], values)
        worksheet.write(row, 3, key[2:], values)
        worksheet.write(row, 4, value, values)

    row: int = 0
    nucleotideList: List[str] = []
    tempNucleotideList: List[str] = []
    for key, value in positionDictionary.items():
        row += 1
        for individualPosition in value:
            tempNucleotideList.append(dnaTable[int(individualPosition)])
            nucleotideString = ','.join(str(codon) for codon in tempNucleotideList)

        nucleotideList.append(nucleotideString)
        worksheet.write(row, 2, nucleotideString, values)
        tempNucleotideList = []

        value = ','.join(str(number) for number in value)
        worksheet.write(row, 5, value, values)

    #closing the file
    workbook.close()

    return nucleotideList

def plotGraphForMutationDatasetAnalysis(fileName, aminoAcidList, mutatedAminoAcidList, nucleotideAnalysisList, mutatedNucleotideAnalysisList, occurencesList, positionList) -> None:

    global RootPath
    #starts the analysis
    os.chdir(RootPath + f'/output/{fileName}/')
    plt.style.use('fivethirtyeight')
    #AMINO ACID
    aminoAcidAnalysisDataFrame = pd.DataFrame({'Amino Acid': aminoAcidList, 'Mutated Amino Acid': mutatedAminoAcidList, 'Occurences': occurencesList, 'Positions': positionList})
    aminoAcidAnalysisDataFrame.plot(title='Amino Acid Dataset Analysis', grid=True, fontsize=10, figsize=(20,6), x=['Amino Acid','Mutated Amino Acid'], y='Occurences', kind='bar', width=0.9, alpha=0.75)
    plt.savefig(f'{fileName}AminoAcid.pdf',bbox_inches='tight')
    #NUCLEOTIDE
    nucleotideAnalysisDataFrame = pd.DataFrame({'Nucleotide': nucleotideAnalysisList, 'Mutated Nucleotide': mutatedNucleotideAnalysisList, 'Occurences': occurencesList, 'Positions': positionList})
    nucleotideAnalysisDataFrame.plot(title='Nucleotide Dataset Analysis', grid=True, fontsize=10, figsize=(20,6), x=['Nucleotide','Mutated Nucleotide'], y='Occurences', kind='bar', width=0.9, alpha=0.75)
    #saves the graph in a file
    plt.savefig(f'{fileName}Nucleotide.pdf',bbox_inches='tight')
