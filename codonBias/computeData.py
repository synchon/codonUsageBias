"""
computes the required data from codons obtained
2016
Synchon Mandal
"""
def codonCompute(codonList, codonTable):

   aminoAcidList = []
   oneLetterList = []

   for codon in codonList:
      for knownCodonKey, knownCodonValue in codonTable.items():
         if codon == knownCodonKey:
            aminoAcidList.append(codon)
            oneLetterList.append(knownCodonValue)

   #removes duplicates
   aminoAcidSet = set(aminoAcidList)
   oneLetterSet = set(oneLetterList)

   setDict = {}
   codeSetDict = {}

   #codon set
   for setValue in aminoAcidSet:
      aminoAcidDuplicates = aminoAcidList.count(setValue)
      setDict[setValue] = aminoAcidDuplicates

   setDict = sorted(setDict.items(), key = lambda t: t[1], reverse = True)

   #amino acid set
   for codeSetValue in oneLetterSet:
      oneLetterDuplicates = oneLetterList.count(codeSetValue)
      codeSetDict[codeSetValue] = oneLetterDuplicates

   codeSetDict = sorted(codeSetDict.items(), key = lambda t: t[1], reverse = True)


   returnList = []
   returnList.append(setDict)
   returnList.append(codeSetDict)
   return returnList
