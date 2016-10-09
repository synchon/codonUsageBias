"""
codonUsageBias
2016
Synchon Mandal
"""
#!/usr/bin/python3

import sys, time, jsonImport, dnaSort, computeData, htmlRender, webbrowser

def main():
    #imports the DNA codon table
    print('Importing codon table...')
    codonTab = jsonImport.openJson('dnaCodonTable')
    print('Codon table imported.\n')

    #runs the sorting function for conversion to triplets
    fileName = str(input('Enter the file name:\t'))
    print('Sorting the genome file...')
    startTimeOfSorting = time.time()
    sortedList = dnaSort.fileSort(fileName)
    totalTimeOfSorting = round(time.time() - startTimeOfSorting, 2)
    print('Sorted.\n')
    print('Total time required for sorting: %f seconds.\n' %(totalTimeOfSorting))

    #runs the computation function for the data obtained
    print('Translating into amino acids...')
    startTimeTranslating = time.time()
    setDict, codeSetDict = computeData.codonCompute(sortedList, codonTab)
    totalTimeTranslating = round(time.time() - startTimeTranslating, 2)
    print('Translated.\n')
    print('Total time required for translating: %f seconds.\n' %(totalTimeTranslating))

    #converts the data into a table in HTML
    htmlPageName = input('Enter title for your output file:\t')
    htmlRender.generateHTMLTable(htmlPageName, setDict, codeSetDict)

    #opens the web browser for viewing the file
    query = input('Press return to open your output file.')
    fileName = 'file:///Users/Synchon/Desktop/codonUsageBias/codonBias' + htmlPageName + '.html'
    webbrowser.open_new(fileName)


if __name__ == '__main__':
    main()
