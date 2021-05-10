import pandas as pd

import sys
import os

script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append("/".join(script_dir.split("/")[:-1]))

from app import db
from app.api.models import HLAAlleles, Peptides, PeptidesPositions, HLAAllelesPeptides


HLA_class = sys.argv[1]  # HLA-I or HLA-II
alleles_list_path = sys.argv[2]
peptides_table_path = sys.argv[3]
affinity_table_path = sys.argv[4]


def fix_HLA_II_name(name):
    if "DRB1_" in name:
        return "HLA-DRB1*" + name[5:7] + ":" + name[7:]
    elif "HLA-DPA1" in name:
        if len(name) == 21:
            return "HLA-DPA1*" + name[8:10] + ":" + name[10:12] + "/DPB1*" + name[17:19] + ":" + name[19:]
        else:
            return "HLA-DPA1*" + name[8:10] + ":" + name[10:12] + "/DPB1*" + name[17:20] + ":" + name[20:]
    else:
        return "HLA-DQA1*" + name[8:10] + ":" + name[10:12] + "/DQB1*" + name[17:19] + ":" + name[19:]

# First, HLA alleles

df = pd.read_csv(alleles_list_path, index_col=0, header=None)
if HLA_class == "HLA-I":
    df.index = df.index.str.replace("HLA-A", "HLA-A*")
    df.index = df.index.str.replace("HLA-B", "HLA-B*")
    df.index = df.index.str.replace("HLA-C", "HLA-C*")
elif HLA_class == "HLA-II":
    df.index = [fix_HLA_II_name(i) for i in df.index]

objects = [
    HLAAlleles(
        hla_allele=hla_allele,
        gene=("DRB1" in hla_allele and "HLA-DRB1") or ("DPA1" in hla_allele and "HLA-DPA1/DPB1") or "HLA-DQA1/DQB1",
        HLA_class=HLA_class.split("-")[1]
    )
    for hla_allele in df.index
]
db.session.bulk_save_objects(objects)
db.session.commit()


# Next, peptides
df = pd.read_csv(peptides_table_path, index_col=0)

objects = [
    Peptides(peptide=peptide)
    for peptide in df.index.unique()
]
db.session.bulk_save_objects(objects)
db.session.commit()

objects = [
    PeptidesPositions(
        peptide=peptide,
        protein=protein,
        start=start,
        end=end,
        gisaid_id="EPI_ISL_402125"
    )
    for peptide, (protein, start, end) in df.iterrows()
]
db.session.bulk_save_objects(objects)
db.session.commit()


# Finally, HLA-peptide affinities
f = open(affinity_table_path)
f.readline()

for line in f:
    hla_allele, peptide, affinity = line.strip().split(",")
    affinity = int(affinity)

    try:
        new_entry = HLAAllelesPeptides(
            hla_allele=fix_HLA_II_name(hla_allele),
            peptide=peptide,
            affinity=min(affinity, 32767)
        )
        db.session.add(new_entry)
        db.session.commit()
        del new_entry
    except Exception as err:
        print(err)
        db.session.rollback()
