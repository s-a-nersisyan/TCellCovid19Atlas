import pandas as pd

import sys
import os

script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append("/".join(script_dir.split("/")[:-1]))

from app import db
from app.api.models import HLAAlleles, Peptides, PeptidesPositions, HLAAllelesPeptides

db.drop_all()
db.create_all()

HLA_class = sys.argv[1]  # HLA-I or HLA-II
alleles_list_path = sys.argv[2]
peptides_table_path = sys.argv[3]
affinity_table_path = sys.argv[4]

# First, HLA alleles
df = pd.read_csv(alleles_list_path, index_col=0, header=None)
if HLA_class == "HLA-I":
    df.index = df.index.str.replace("HLA-A", "HLA-A*")
    df.index = df.index.str.replace("HLA-B", "HLA-B*")
    df.index = df.index.str.replace("HLA-C", "HLA-C*")

objects = [
    HLAAlleles(
        hla_allele=hla_allele,
        gene=hla_allele[0:5],
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
df = pd.read_csv(affinity_table_path)
df = df.drop_duplicates()

for _, (hla_allele, peptide, affinity) in df.iterrows():
    new_entry = HLAAllelesPeptides(
        hla_allele=hla_allele,
        peptide=peptide,
        affinity=min(affinity, 32767)
    )
    db.session.add(new_entry)
    db.session.commit()
    del new_entry

'''

objects = [
    HLAAllelesPeptides(
        hla_allele=hla_allele,
        peptide=peptide,
        affinity=affinity
    )
    for _, (hla_allele, peptide, affinity) in df.iterrows()
]
print("Objects generated!")
db.session.bulk_save_objects(objects)
db.session.commit()
'''
