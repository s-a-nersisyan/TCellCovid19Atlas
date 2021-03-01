import pandas as pd
import tqdm
import sys
import os

peptide = sys.argv[1]
affinity_file = open(peptide + ".aff", "r")
tb = pd.DataFrame(columns=["HLA", "affinity"])
iterator = iter(affinity_file)
for line in iterator:
    if not ("Pos" in line.split(" ")):
        continue
    line = next(iterator)
    line = next(iterator)
    print(line.split()[1], line.split()[14])
    tb.loc[len(tb)] = [
        line.split()[1],
        line.split()[14]
    ]

tb.to_csv("../affinities/" + peptide + ".tsv", sep="\t", index=None)
affinity_file.close()
