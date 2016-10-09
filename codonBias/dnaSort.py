"""
sorts the DNA and returns the sorted file
2016
Synchon Mandal
"""
import math

#this is currently for FASTA format only
def fileSort(fileName):

   with open('%s.txt' %(str(fileName)),'rU') as genomeFile:
      next(genomeFile)#skips the first line of FASTA file
      sortedFile = []
      readGenome = genomeFile.read()
      readGenome = str(readGenome)
      readGenome = readGenome.split('\n')
      readGenome = ''.join(readGenome)
      iterVar = 0

      if len(readGenome)%3 == 0:
         for codon in range(int(len(readGenome)/3)):
            triplet = readGenome[iterVar] + readGenome[iterVar + 1] + readGenome[iterVar + 2]
            sortedFile.append(triplet)
            iterVar += 3

      elif len(readGenome)%3 != 0:
         for codon in range(int(math.floor(len(readGenome)/3))):
            try:
                triplet = readGenome[iterVar] + readGenome[iterVar + 1] + readGenome[iterVar + 2]
                sortedFile.append(triplet)
                iterVar += 3
            except:
                pass

   print ('Total number of codons: %i' %(len(sortedFile)))
   return sortedFile
