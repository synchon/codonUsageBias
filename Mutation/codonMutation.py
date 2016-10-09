"""
mutation using codon variables
2016
Synchon Mandal
"""
import jsonImport

def main():
    codon = str(input("Enter a codon:\n"))
    aminoAcidCodon, aminoAcidBase, tripletVar = codonIdentVar(list(codon))
    print('%s --> %s[%s]' %(aminoAcidCodon, aminoAcidBase, tripletVar))
    print("\nMutated list:\n")
    codonMutation(codon)

def codonIdentVar(codonList):
   #import codon table
   codonTable = jsonImport.openJson('dnaCodonTable')

   #define bases
   bases = {'A': 'alpha1','G': 'alpha2','T':'beta1','C': 'beta2'}


   #join to form a string
   #param 1
   codonString = ''.join(codonList)

   baseVar = []
   #param 2
   aminoAcid = ''

   #run loop for identification
   for codon in codonList:
       for baseKey, baseValue in bases.items():
           if codon == baseKey:
               baseVar.append(baseValue)
   #param 3
   baseVarString = ','.join(baseVar)

   #run the check for amino acid identification
   for knownCodonKey, knownCodonValue in codonTable.items():
       if codonString == knownCodonKey:
           aminoAcid = str(knownCodonValue)

   returnList = []
   returnList.append(codonString)
   returnList.append(aminoAcid)
   returnList.append(baseVarString)
   return returnList

def codonMutation(codon):
     #define bases
     bases = {'A': 'alpha1','G': 'alpha2','T':'beta1','C': 'beta2'}

     codonList = list(codon)
     #starts at first and iterates over
     for index, value in enumerate(codonList):
         for baseKey, baseValue in bases.items():
             if baseKey != value:
                codonList[index] = baseKey
                codonIdentVar(codonList)
                aminoAcidCodon, aminoAcidBase, tripletVar = codonIdentVar(codonList)
                print('%s --> %s[%s]' %(aminoAcidCodon, aminoAcidBase, tripletVar))
                print('=======================================')
                codonList = list(codon)


if __name__ == '__main__':
     main()
