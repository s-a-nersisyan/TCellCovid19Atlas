import pandas as pd
import tqdm
import sys
import os


# Path to FASTA file with viral proteome
proteome_path = sys.argv[1]
# Comma-separated list of peptide lengths
peptide_lengths = list(map(int, sys.argv[2].split(",")))
# Output: comma-separated table of the following form:
# Peptide, protein name, start coordinate, end coordinate
# Coordinate system is 0-based

proteome_file = open(proteome_path)
peptides = set()
for line in proteome_file:
    if line.startswith(">"):
        protein_name = line.split("|")[0].lstrip(">")
    else:
        protein = line.strip().rstrip("*")
        for k in peptide_lengths:
            for start in range(len(protein) - k + 1):
                peptides.add(protein[start:start+k])

print("Peptide,Protein,Start,End")
for peptide in peptides:
    print("{},{},{},{}".format(
        peptide,
        "-",
        "-",
        "-"
    ))
