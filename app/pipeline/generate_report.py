import pandas as pd
import numpy as np
from natsort import natsorted

import json
import sys
import os

script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append("/".join(script_dir.split("/")[:-2]))

from app import db
from app.api.models import HLAAllelesPeptides

gisaid_id = sys.argv[1]

AFFINITY_THRESHOLD1 = 500
AFFINITY_THRESHOLD2 = 500
peptide_len_range = {
    "HLA-I": [8, 14],
    "HLA-II": [15, 20]
}

df = pd.read_csv("{}/diff_HLA-I.csv".format(gisaid_id))
#df = pd.read_csv("{}/binding_affinities_HLA-I.csv".format(gisaid_id))
#print(df)

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
    
    block_started = False
    for i, (r, m) in enumerate(zip(ref_seq, mut_seq)):
        if r == m:
            if block_started:
                block_started = False
                block_end = i
                mutations.append([ref_protein, block_start, block_end])
        else:
            if not block_started:
                block_started = True
                block_start = i
    
    if block_started:
        block_started = False
        block_end = i
        mutations.append([ref_protein, block_start, block_end])

report = []
for protein, start, end in mutations:
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

        try:
            mut_hla_pep_iter = HLAAllelesPeptides.query.filter_by(peptide=mut_pep).order_by(HLAAllelesPeptides.hla_allele).all()
        except:
            mut_hla_pep_iter = iter(lambda: None, 1)
            db.session.rollback()

        for ref_hla_pep, mut_hla_pep in zip(ref_hla_pep_iter, mut_hla_pep_iter):
            hla_allele = ref_hla_pep.hla_allele if ref_hla_pep else None
            if not hla_allele:
                hla_allele = mut_hla_pep.hla_allele

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
    report.append([protein, start, end, report_row])

report_file = open("{}/report.json".format(gisaid_id), "w")
json.dump(report, report_file)
report_file.close()

df = pd.DataFrame(columns=[
    "Protein", "Aln start", "Aln end", "Allele", "Ref peptide", "Mut peptide", "Ref aff", "Mut aff"
])
for mutation in report:
    df1 = pd.DataFrame(mutation[3])
    df1.columns = ["Allele", "Ref peptide", "Mut peptide", "Ref aff", "Mut aff"]
    df1["Protein"] = mutation[0]
    df1["Aln start"] = mutation[1]
    df1["Aln end"] = mutation[2]
    df = pd.concat([df, df1])

df.to_csv("{}/report.csv".format(gisaid_id), index=None)
