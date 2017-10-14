#!/usr/bin/sh

# Initiate directories:
if [ ! -d "input" ]; then
  mkdir input
fi
if [ ! -d "output" ]; then
  mkdir output
fi

# Download dependencies:
conda install biopython matplotlib pandas
