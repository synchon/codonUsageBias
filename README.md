# Codon Usage Bias and Mutation Prediction
---
## Requires python 3.x

1. Installing package:
	Open a terminal and type `sh install.sh`

2. Initiating input files:
    1. For online use:
	      In the /input directory add *gene*.csv (mutation dataset)

    2. For offline use:
        In the /input directory add *gene*.fasta (gene file) and *gene*.csv (mutation dataset)

3. Running the program:
	./main.py -a 12AB.1 -e synchon@email.com

4. Collecting output:
	Check the /output for the required directory.

NOTE:
### The program automatically determines the required input source i.e., if .fasta file is present locally, it uses that or else searches the online Entrez database.
### In this version, it has no option for being verbose; will update.