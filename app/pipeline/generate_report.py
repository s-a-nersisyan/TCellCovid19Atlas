import pandas as pd
import numpy as np
from natsort import natsorted
from plot_mutations import protein_plot, allele_plot

import json
import sys
import os

script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append("/".join(script_dir.split("/")[:-2]))

from app import db
from app.api.models import HLAAllelesPeptides, HLAAlleles

PROTEINS = [
    "Spike", "N", "M", "E", "NS3", "NS6",
    "NS7a", "NS7b", "NS8", "NS9b", "NS9c",
    "NSP1", "NSP2", "NSP3", "NSP4", "NSP5",
    "NSP6", "NSP7", "NSP8", "NSP9", "NSP10",
    "NSP11", "NSP12", "NSP13", "NSP14", "NSP15", "NSP16"
]

PROTEIN_SIGNIFICANCE = dict(zip(
    PROTEINS,
    range(len(PROTEINS))
))


def get_mutation(ref_seq, mut_seq, block_start, block_end, ref_numspaces, sub_type):
    if (sub_type == "sub"):
        if (block_end - block_start == 1):
            mut = "{}{}{}".format(
                ref_seq[block_start],
                block_start - ref_numspaces + 1,
                mut_seq[block_start]
            )
        else:
            mut = "{}-{} {}->{}".format(
                block_start - ref_numspaces + 1,
                block_end - ref_numspaces,
                ref_seq[block_start:block_end],
                mut_seq[block_start:block_end]
            )
    
    if (sub_type == "del"):
        if (block_end - block_start == 1):
            mut = "{}{} {}".format(
                ref_seq[block_start],
                block_start - ref_numspaces + 1,
                "deletion"
            )
        elif (block_end - block_start <= 5):
            mut = "{} {}-{} {}".format(
                ref_seq[block_start:block_end],
                block_start - ref_numspaces + 1,
                block_end - ref_numspaces,
                "deletion"
            )
        elif (block_end == len(ref_seq) and block_end - block_start > 5):
            mut = "{}{}{}".format(
                ref_seq[block_start],
                block_start - ref_numspaces + 1,
                "stop"
            )
        else:
            mut = "{}-{} {}".format(
                block_start - ref_numspaces + 1,
                block_end - ref_numspaces,
                "deletion"
            )

    if (sub_type == "ins"):
        if (block_end - block_start == 1):
            mut = "{}{} {}".format(
                mut_seq[block_start],
                block_start - ref_numspaces + 1,
                "insertion"
            )
        elif (block_end - block_start <= 5):
            mut = "{} {}-{} {}".format(
                mut_seq[block_start:block_end],
                block_start - ref_numspaces + 1,
                block_end - ref_numspaces,
                "insertion"
            )
 
        elif (block_start == 0 and block_end - block_start > 5):
            mut = "{}{}{}".format(
                mut_seq[block_start],
                block_start - ref_numspaces + 1,
                "start"
            )
        else:
            mut = "{}-{} {}".format(
                block_start - ref_numspaces + 1,
                block_end - ref_numspaces,
                "insertion"
            )
    
    return mut

def get_sub_type(r, m):
    if (r == "-"):
        return "ins"
    if (m == "-"):
        return "del"
    if (r != m):
        return "sub"
    return "idn"

gisaid_id = sys.argv[1]

AFFINITY_THRESHOLD1 = 500
AFFINITY_THRESHOLD2 = 50

peptide_len_range = {
    "HLA-I": [8, 14],
    "HLA-II": [15, 20]
}

# The following block was added by dude
for hla_class in ["HLA-I", "HLA-II"]:
    df = pd.read_csv("{}/diff_{}.csv".format(gisaid_id, hla_class))

    ref_aln = open("{}/ref_aln.fasta".format(gisaid_id))
    mut_aln = open("{}/mut_aln.fasta".format(gisaid_id))

    mutations = []
    alignments = {}
    while 1:
        ref_protein = ref_aln.readline().strip().lstrip(">")
        mut_protein = mut_aln.readline().strip().lstrip(">")
        if not ref_protein and not mut_protein:
            break

        ref_seq = ref_aln.readline().strip()
        mut_seq = mut_aln.readline().strip()

        alignments[ref_protein] = (ref_seq, mut_seq)
        ref_numspaces = 0

        i = 0
        while (i < len(ref_seq)):
            r = ref_seq[i]
            m = mut_seq[i]
            
            if (r == "-"):
                ref_numspaces += 1
     
            if (r == m):
                i += 1
                continue
            
            block_start = block_end = i
            sub_type = get_sub_type(r, m)
            while (block_end < len(ref_seq)):
                r = ref_seq[block_end]
                m = mut_seq[block_end]
                
                if block_end > block_start and r == "-":
                    ref_numspaces += 1
                
                cur_type = get_sub_type(r, m)

                if (cur_type != sub_type):
                    break
                
                block_end += 1
            
            mutation = get_mutation(ref_seq, mut_seq, block_start,
                    block_end, ref_numspaces, sub_type)
            
            ref_block = ref_seq[max(0, block_start - 25):min(len(ref_seq), block_end + 25)]
            mut_block = mut_seq[max(0, block_start - 25):min(len(mut_seq), block_end + 25)]

            mutations.append([ref_protein, block_start, block_end, mutation, ref_block, mut_block])

            i = block_end
            
    report = []
    for protein, start, end, mut, ref_block, mut_block in mutations:
        df1 = df.loc[
            ((df["protein"] == protein) & (df["aln_start"] < end) & (df["aln_end"] > start))
        ]
        df1 = df1.sort_values(["aln_start", "aln_end"])
        report_row = []
        for i, row in df1.iterrows():
            ref_pep = row["ref_pep"]
            mut_pep = row["mut_pep"]
            
            try:
                ref_hla_pep_iter = HLAAllelesPeptides.query.filter_by(peptide=ref_pep).order_by(HLAAllelesPeptides.hla_allele).all()
            except:
                ref_hla_pep_iter = iter(lambda: None, 1)  # Infinite iterator
                db.session.rollback()
            
            if (
                pd.isna(ref_pep) or
                len(ref_hla_pep_iter) == 0 or
                len(ref_pep) < peptide_len_range[hla_class][0] or
                len(ref_pep) > peptide_len_range[hla_class][1]
            ):
                ref_hla_pep_iter = iter(lambda: None, 1)  # Infinite iterator

            try:
                mut_hla_pep_iter = HLAAllelesPeptides.query.filter_by(peptide=mut_pep).order_by(HLAAllelesPeptides.hla_allele).all()
            except:
                mut_hla_pep_iter = iter(lambda: None, 1)
                db.session.rollback()
            
            if (
                pd.isna(mut_pep) or
                len(mut_hla_pep_iter) == 0 or
                len(mut_pep) < peptide_len_range[hla_class][0] or
                len(mut_pep) > peptide_len_range[hla_class][1]
            ):
                mut_hla_pep_iter = iter(lambda: None, 1)  # Infinite iterator

            for ref_hla_pep, mut_hla_pep in zip(ref_hla_pep_iter, mut_hla_pep_iter):
                if (ref_hla_pep != None and mut_hla_pep != None):
                    hla_allele = ref_hla_pep.hla_allele
                    
                    try:
                        hla_flag = HLAAlleles.query.filter(
                            HLAAlleles.hla_allele == ref_hla_pep.hla_allele
                        ).filter(HLAAlleles.HLA_class == hla_class.split("-")[-1]).count()
                            
                        hla_flag = hla_flag and HLAAlleles.query.filter(
                            HLAAlleles.hla_allele == mut_hla_pep.hla_allele
                        ).filter(HLAAlleles.HLA_class == hla_class.split("-")[-1]).count()
                    except:
                        hla_flag = 0
                        db.session.rollback()
                else:
                    hla_allele = ref_hla_pep.hla_allele if ref_hla_pep else None
                    if not hla_allele:
                        hla_allele = mut_hla_pep.hla_allele

                    try:
                        hla_flag = HLAAlleles.query.filter(
                            HLAAlleles.hla_allele == hla_allele
                        ).filter(HLAAlleles.HLA_class == hla_class.split("-")[-1]).count()
                    except:
                        hla_flag = 0
                        db.session.rollback()

                if not hla_flag:
                    print(hla_allele, hla_class.split("-")[-1])
                    continue

                ref_aff = ref_hla_pep.affinity if ref_hla_pep else None
                mut_aff = mut_hla_pep.affinity if mut_hla_pep else None

                if not ref_aff and mut_aff > AFFINITY_THRESHOLD1:
                    continue
                if not mut_aff and ref_aff > AFFINITY_THRESHOLD1:
                    continue

                #if ref_aff and mut_aff and \
                #    not ((ref_aff <= AFFINITY_THRESHOLD1 and mut_aff > AFFINITY_THRESHOLD2) or \
                #    (ref_aff > AFFINITY_THRESHOLD2 and mut_aff <= AFFINITY_THRESHOLD1)):
                #    continue
                
                if ref_aff and mut_aff and ref_aff > AFFINITY_THRESHOLD1 and mut_aff > AFFINITY_THRESHOLD1:
                    continue

                if pd.isna(ref_pep):
                    ref_pep = None
                if pd.isna(mut_pep):
                    mut_pep = None

                if ref_aff and mut_aff and np.abs(np.log2(ref_aff/mut_aff)) < 1:
                    continue
                
                report_row.append([
                    hla_allele, ref_pep, mut_pep, ref_aff, mut_aff
                ])
        
        def sorting_key(row):
            if row[3] and row[4]:
                neg_ratio = -np.abs(np.log2(row[3]/row[4]))
                min_aff = min(row[3], row[4])
            else:
                neg_ratio = -np.inf
                if row[3]:
                    min_aff = row[3]
                else:
                    min_aff = row[4]

            return neg_ratio, min_aff, row[0]
            
        report_row = sorted(report_row, key=sorting_key)
        report.append([protein, start, end, mut, report_row])
    
    report_file = open("{}/report_{}.json".format(gisaid_id, hla_class), "w")
    json.dump(report, report_file)
    report_file.close()
    
    # report_file = open("{}/report_{}.json".format(gisaid_id, hla_class), "r")
    # report = json.load(report_file)
    # report_file.close()

    df = pd.DataFrame(columns=[
        "Protein", "Aln start", "Aln end",
        "Mutation", "Allele", "Ref peptide",
        "Mut peptide", "Ref aff", "Mut aff",
        "Ref aff type", "Mut aff type"
    ])

    for mutation in report:
        df1 = pd.DataFrame(mutation[4])
        if df1.empty:
            continue
        else:
            print(df1)
        df1.columns = ["Allele", "Ref peptide", "Mut peptide", "Ref aff", "Mut aff"]
        df1["Protein"] = mutation[0]
        df1["Aln start"] = mutation[1]
        df1["Aln end"] = mutation[2]
        df1["Mutation"] = mutation[3]
        df1["Ref aff type"] = df1["Mut aff type"] = None
        
        for i, row in df1.iterrows(): 
            if (pd.isna(row["Ref aff"]) or
                    row["Ref aff"] >= AFFINITY_THRESHOLD1):
                df1.loc[i, "Ref aff type"] = "Weak/no binding"
            elif (row["Ref aff"] >  AFFINITY_THRESHOLD2 and
                    row["Ref aff"] <= AFFINITY_THRESHOLD1):
                df1.loc[i, "Ref aff type"] = "Moderate binding" 
            elif (row["Ref aff"] <= AFFINITY_THRESHOLD2):
                df1.loc[i, "Ref aff type"] = "Tight binding"

            if (pd.isna(row["Mut aff"]) or
                    row["Mut aff"] >= AFFINITY_THRESHOLD1):
                df1.loc[i, "Mut aff type"] = "Weak/no binding"
            elif (row["Mut aff"] >  AFFINITY_THRESHOLD2 and
                    row["Mut aff"] <= AFFINITY_THRESHOLD1):
                df1.loc[i, "Mut aff type"] = "Moderate binding" 
            elif (row["Mut aff"] <= AFFINITY_THRESHOLD2):
                df1.loc[i, "Mut aff type"] = "Tight binding"

        df = pd.concat([df, df1])
    
    df.to_csv("{}/report_{}.csv".format(gisaid_id, hla_class), index=None)

    # create tables dir    
    if not os.path.exists("../static/tables/{}/{}/".format(gisaid_id, hla_class)):
        os.makedirs("../static/tables/{}/{}/".format(gisaid_id, hla_class))
    
    # create plots dir    
    if not os.path.exists("../static/plots/{}/{}/".format(gisaid_id, hla_class)):
        os.makedirs("../static/plots/{}/{}/".format(gisaid_id, hla_class))

    # generate allele plots
    alleles = set(df["Allele"])

    for allele in sorted(alleles):
        print(gisaid_id, allele, hla_class)
        allele_plot(gisaid_id, hla_class, allele)

    # generate protein plots
    proteins = set(df["Protein"])
    proteins = sorted(proteins, key=lambda protein: PROTEIN_SIGNIFICANCE[protein])
    for protein in proteins + ["Summary"]:
        print(gisaid_id, protein, hla_class)
        protein_plot(gisaid_id, hla_class, protein)

    # generate mutation summary table
    if (hla_class == "HLA_II"):
        continue
    
    proteins = pd.DataFrame(proteins, columns=["Protein"])
    proteins.to_csv("{}/proteins.csv".format(gisaid_id), index=None)

    mut_df = pd.DataFrame(columns=[
        "Protein", "Block start",
        "Block end", "Mutation",
        "Ref block", "Mut block"
    ])

    proteins = set(df["Protein"])
    for protein, start, end, mut, ref_block, mut_block in mutations:
        if not protein in proteins:
            continue

        mut_df.loc[len(mut_df)] = [
            protein,
            start - (max(0, start - 25)),
            end - (max(0, start - 25)),
            mut, ref_block, mut_block
        ]

    mut_df.to_csv("{}/mutations.csv".format(gisaid_id), index=None)
