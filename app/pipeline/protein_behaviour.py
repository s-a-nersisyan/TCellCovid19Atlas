import pandas as pd
import numpy as np

import sys
import os

import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams["font.sans-serif"] = "Arial"
mpl.rcParams["figure.dpi"] = 300

def get_mutation_identifier(gisaid_id, protein):
    if (protein == None):
        return None

    ref_aln = None
    ref_aln_file = open("{}/ref_aln.fasta".format(gisaid_id), "r")
    ref_iter = iter(ref_aln_file)
    for line in ref_iter:
        if (line.startswith(">" + protein)):
            ref_aln = next(ref_iter)
            break
    if (ref_aln == None):
        return None
    ref_aln_file.close()

    mut_aln = None
    mut_aln_file = open("{}/mut_aln.fasta".format(gisaid_id), "r")
    mut_iter = iter(mut_aln_file)
    for line in mut_aln_file:
        if (line.startswith(">" + protein)):
            mut_aln = next(mut_iter)
            break
    if (mut_aln == None):
        return None
    mut_aln_file.close()

    result = protein
    for ind in range(len(ref_aln)):
        if (ref_aln[ind] == "-"):
            gap_num += 1

        if (ref_aln[ind] != mut_aln[ind]):
            result += ref_aln[ind] + str(ind) + mut_aln[ind]

    return result


gisaid_id = sys.argv[1]
if (len(sys.argv) > 2):
    protein = sys.argv[2]
else:
    protein = None

report = pd.read_csv("{}/report.csv".format(
    gisaid_id
))

if (protein != None):
    report = report[report["Protein"] == protein]
report = report.fillna(np.inf)

alleles = set(report["Allele"])

df = pd.DataFrame(columns=["Allele", "Weaker binding", "Stronger binding"])
for allele in sorted(alleles):
    increase = len(report[
        (report["Allele"] == allele) &
        (report["Ref aff"] < report["Mut aff"])
    ])

    decrease = len(report[
        (report["Allele"] == allele) &
        (report["Ref aff"] > report["Mut aff"])
    ])

    df.loc[len(df)] = [allele, increase, decrease]

df = df.set_index(df["Allele"])
df = df.sort_values(by=["Weaker binding", "Stronger binding"], ascending=True)
max_val = np.max(df["Weaker binding"])
max_val = max(
    np.max(df["Weaker binding"]),
    np.max(df["Stronger binding"])
) + 1
df["Stronger binding"] = -df["Stronger binding"]

sns.set()
df.plot.barh(
    stacked=True,
    xlabel="", ylabel="",
    title=protein,
    figsize=(6, len(df) / 4)
)
plt.xlim([-max_val, max_val])
plt.tight_layout()
plt.legend(loc='lower right')

plt.savefig(
    os.path.join("{}/plots/{}.png".format(
        gisaid_id,
        protein
    ))
)
plt.close()





