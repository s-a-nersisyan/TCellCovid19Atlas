import numpy as np
import pandas as pd
import sys
import os

# alignment tool
from Bio import pairwise2
from Bio.pairwise2 import format_alignment

MATCH_COST = 5
MISMATCH_COST = -4
GAP_COST = -10
GAP_EXTEND_COST = -0.5

ref_fasta_path = sys.argv[1]
mut_fasta_path = sys.argv[2]
ref_pep_path = sys.argv[3]
mut_pep_path = sys.argv[4]
output_dir = sys.argv[5]

def get_epi(path):
    start = path.rfind("EPI_ISL_") + len("EPI_ISL_")
    end = start + path[start:].find("/")
    return path[start:end]

def protein_dict(fasta_path):
    protein_dict = dict()
    fasta_file = open(fasta_path, "r")
    iterator = iter(fasta_file)

    for line in iterator:
        if not (line.startswith(">")):
            continue

        line = line.lstrip(">")
        protein_name = line.split("|")[0]
        protein = next(iterator)

        protein_dict[protein_name] = protein

    fasta_file.close()
    return protein_dict

def protein_alignment_dict(ref_protein_dict, mut_protein_dict):
    protein_alignment_dict = dict()
    for protein_name in ref_protein_dict:
        ref_protein = ref_protein_dict[protein_name]
        mut_protein = mut_protein_dict[protein_name]

        alignment = [*pairwise2.align.globalms(
            ref_protein, mut_protein,
            MATCH_COST, MISMATCH_COST,
            GAP_COST, GAP_EXTEND_COST
        )[0]]

        protein_alignment_dict[protein_name] = (
            alignment[0].rstrip("*\n"),
            alignment[1].rstrip("*\n")
        )

    return protein_alignment_dict

def conversion(reference, alignment, direction=None):
    forward = list()
    backward = list()

    gap_count = 0
    for ind, nucl in enumerate(alignment):
        if (nucl != "-"):
            forward.append(ind)
        else:
            gap_count += 1
        backward.append(ind - gap_count)

    if (direction == "forward"):
        return forward
    if (direction == "backward"):
        return backward
    return forward, backward

ref_protein_dict = protein_dict(ref_fasta_path)
mut_protein_dict = protein_dict(mut_fasta_path)
protein_names = [protein_name for protein_name in ref_protein_dict]

# check: all referece proteins are included in mutation fasta file
for protein_name in protein_names:
    if not (protein_name in mut_protein_dict):
        print("Err (diff): protein", protein_name)

alignments = protein_alignment_dict(
    ref_protein_dict,
    mut_protein_dict
)

ref_pep_tb = pd.read_csv(ref_pep_path, sep=",")
mut_pep_tb = pd.read_csv(mut_pep_path, sep=",")

output_tb = pd.DataFrame(columns=[
    "protein",
    "ref_pep",
    "mut_pep",
    "ref_aln_pep",
    "mut_aln_pep",
    "ref_ind",
    "mut_ind",
    "aln_ind"
])

for protein_name in protein_names:
    print(protein_name)

    ref_tb = ref_pep_tb[
        ref_pep_tb["Protein"] == protein_name
    ]

    mut_tb = mut_pep_tb[
        mut_pep_tb["Protein"] == protein_name
    ]

    diff_pep = set(ref_tb["Peptide"]) - set(mut_tb["Peptide"])

    ref_tb = ref_tb.set_index("Peptide")
    mut_tb = mut_tb.set_index("Peptide")

    for peptide in diff_pep:
        start = ref_tb.loc[peptide]["Start"]
        end = ref_tb.loc[peptide]["End"]

        ref_align, mut_align = alignments[protein_name]
        fw = conversion(ref_protein_dict[protein_name], ref_align,
                        direction="forward")

        mut_protein = mut_protein_dict[protein_name]
        bw = conversion(mut_protein, mut_align,
                        direction="backward")

        output_tb.loc[len(output_tb)] = [
            protein_name,
            peptide,
            mut_protein[bw[fw[start]]:bw[fw[end - 1]] + 1],
            ref_align[fw[start]:fw[end - 1] + 1],
            mut_align[fw[start]:fw[end - 1] + 1],
            "{}:{}".format(start, end),
            "{}:{}".format(bw[fw[start]], bw[fw[end - 1]] + 1),
            "{}:{}".format(fw[start], fw[end - 1] + 1)
        ]

output_path = os.path.join(output_dir, "{}_{}/".format(
    get_epi(ref_fasta_path),
    get_epi(mut_fasta_path)
))

print(output_path)

os.mkdir(output_path)
ref_aln_fasta = open(os.path.join(output_path, "ref_aln.fasta"), "w")
mut_aln_fasta = open(os.path.join(output_path, "mut_aln.fasta"), "w")

for protein_name in alignments:
    ref_aln_fasta.write(">" + protein_name + "\n")
    ref_aln_fasta.write(alignments[protein_name][0] + "\n")
    mut_aln_fasta.write(">" + protein_name + "\n")
    mut_aln_fasta.write(alignments[protein_name][1] + "\n")

output_tb.to_csv(os.path.join(output_path, "diff.csv"), index=None, sep=",")
