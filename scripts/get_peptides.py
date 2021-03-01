import pandas as pd
import tqdm
import sys
import os

PROTEINS_PATH = "/home/steve/huge/steve/HLA/data/GISAID/allprot0918.fasta"
gisaid_id= sys.argv[1].lower()

peptide_lens = [8, 9, 10, 11, 12, 13, 14]

def get_peptides(protein):
    peptides = []
    for start in range(len(protein)):
        for length in peptide_lens:
            if (start + length < len(protein)):
                peptides.append((
                    protein[start:start + length],
                    start,
                    start + length
                ))
    return peptides


tb = pd.DataFrame(columns=["peptide", "protein", "begin", "end"])

gisaid_fasta = open(PROTEINS_PATH, "r")
iterator = iter(gisaid_fasta)

for ind, description in enumerate(tqdm.tqdm(iterator)):
    protein = next(iterator).rstrip("\n")
    line_id = description.split("|")[3].lower()
    name = description.split("|")[0].lstrip(">")

    if not (line_id == gisaid_id):
        continue

    peptides = get_peptides(protein)
    for peptide in peptides:
        tb.loc[len(tb)] = [peptide[0], name,
                           peptide[1], peptide[2]]

    if (ind % 10 == 0):
        tb.to_csv("../" + gisaid_id.upper() + "_peptides.tsv", sep="\t", index=None)

tb.to_csv("../" + gisaid_id.upper() + "_peptides.tsv", sep="\t", index=None)
gisaid_fasta.close()

