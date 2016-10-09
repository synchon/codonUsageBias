"""
codon identification using variables
2016
Synchon
"""
import sys, random, jsonImport

def main():
    print('Generating a random codon:\n')
    aminoAcidCodon, aminoAcidBase, tripletVar = codonIdentVar()
    print('%s --> %s[%s]' %(aminoAcidCodon, aminoAcidBase, tripletVar))


def codonIdentVar():
   #import codon table
   codonTab = jsonImport.openJson('dnaCodonTable')

   #define bases
   bases = {'A': 'alpha1','C': 'alpha2', 'G': 'beta1', 'T':'beta2'}

   #generate a random triplet
   codonList = random.sample(list(bases.keys()),3)

   #join to form a string
   #param 1
   codonStr = ''.join(codonList)

   baseVar = []
   #param 2
   aminoAcid = ''

   #run loop for identification
   for codon in codonList:
       for baseKey, baseValue in bases.items():
           if codon == baseKey:
               baseVar.append(baseValue)
   #param 3
   baseVarStr = ','.join(baseVar)

   #run the check for amino acid identification
   for knownCodonKey, knownCodonValue in codonTab.items():
       if codonStr == knownCodonKey:
           aminoAcid = str(knownCodonValue)

   returnList = []
   returnList.append(codonStr)
   returnList.append(aminoAcid)
   returnList.append(baseVarStr)
   return returnList

if __name__ == "__main__":
    main()
