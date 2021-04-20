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


def conversion(reference, alignment):
    return [ind for ind, nucl in enumerate(alignment) if nucl != "-"]


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
    "ref_start",
    "ref_end",
    "mut_start",
    "mut_end",
    "aln_start",
    "aln_end"
])

for protein_name in protein_names:
    ref_tb = ref_pep_tb[
        ref_pep_tb["Protein"] == protein_name
    ]

    mut_tb = mut_pep_tb[
        mut_pep_tb["Protein"] == protein_name
    ]

    diff_pep = set(ref_tb["Peptide"]) - set(mut_tb["Peptide"])

    ref_tb = ref_tb.set_index("Peptide")
    mut_tb = mut_tb.set_index("Peptide")

    for ref_peptide in diff_pep:
        ref_start = int(ref_tb.loc[ref_peptide]["Start"])
        ref_end = int(ref_tb.loc[ref_peptide]["End"])

        ref_align, mut_align = alignments[protein_name]
        coord = conversion(ref_protein_dict[protein_name], ref_align)
        aln_start, aln_end = coord[ref_start], coord[ref_end - 1] + 1
        
        ref_peptide_aln = ref_align[aln_start:aln_end]
        mut_peptide_aln = mut_align[aln_start:aln_end]

        mut_peptide = mut_peptide_aln.replace("-", "")
        if mut_peptide in ref_tb.index:
            mut_peptide = ""

        if mut_peptide in mut_tb.index:
            mut_start = int(mut_tb.loc[mut_peptide]["Start"])
            mut_end = int(mut_tb.loc[mut_peptide]["End"])
        else:
            mut_start, mut_end = None, None

        output_tb.loc[len(output_tb)] = [
            protein_name,
            ref_peptide, mut_peptide,
            ref_peptide_aln, mut_peptide_aln,
            ref_start, ref_end,
            mut_start, mut_end,
            aln_start, aln_end
        ]

ref_aln_fasta = open(os.path.join(output_dir, "ref_aln.fasta"), "w")
mut_aln_fasta = open(os.path.join(output_dir, "mut_aln.fasta"), "w")

for protein_name in alignments:
    ref_aln_fasta.write(">" + protein_name + "\n")
    ref_aln_fasta.write(alignments[protein_name][0] + "\n")
    mut_aln_fasta.write(">" + protein_name + "\n")
    mut_aln_fasta.write(alignments[protein_name][1] + "\n")

output_tb.to_csv(os.path.join(output_dir, "diff.csv"), index=None, sep=",")
