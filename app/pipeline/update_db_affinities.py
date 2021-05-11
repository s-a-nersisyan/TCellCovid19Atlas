import pandas as pd

import sys
import os

script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append("/".join(script_dir.split("/")[:-2]))

from app import db
from app.api.models import Peptides, HLAAllelesPeptides

gisaid_id = sys.argv[1]

# This function is a dublicate of one
# taken from TCell../db/fill_HLA_peptide_tables.py
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

for hla_class in ["HLA-I", "HLA-II"]:
    df = pd.read_csv("{}/novel_peptides_{}.csv".format(gisaid_id, hla_class))
    
    for peptide in df["Peptide"]:
        try:
            new_entry = Peptides(peptide=peptide)
            db.session.add(new_entry)
            db.session.commit()
            del new_entry
        except Exception as err:
            db.session.rollback()
    
    df = pd.read_csv("{}/binding_affinities_{}.csv".format(gisaid_id, hla_class))
    for _, (hla_allele, peptide, affinity) in df.iterrows():
        if (hla_class == "HLA-II"):
            hla_allele = fix_HLA_II_name(hla_allele)

        try:
            new_entry = HLAAllelesPeptides(
                hla_allele=hla_allele,
                peptide=peptide,
                affinity=min(affinity, 32767)
            )
            db.session.add(new_entry)
            db.session.commit()
            del new_entry
        except:
            db.session.rollback()
