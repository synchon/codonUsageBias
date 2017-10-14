# Codon Usage Bias and Mutation Prediction
### Requires python 3.6 (Conda distribution recommended)

1. Installing Conda package manager and dependencies:
	First, check if you have Conda installed, by typing `conda -v` in your terminal. If you don't have it, please go to <https://www.continuum.io/downloads> and download the Python 3.6 version. After that, type `sh install.sh` to install the dependencies and initiate the program.

2. Initiating input files:
    1. For online use:
	      In the /input directory add `*gene*.csv` (mutation dataset)

    2. For offline use:
        In the /input directory add `*gene*.fasta` (gene file) and `*gene*.csv` (mutation dataset)

3. Running the program:
	python3 main.py -a=CR541742.1 -e=abc@email.com

4. Collecting output:
	Check the /output for the required directory.

NOTE:
#### The program automatically determines the required input source i.e., if .fasta file is present locally, it uses that or else searches the online Entrez database.
#### In this version, it has no option for being verbose; will update.
