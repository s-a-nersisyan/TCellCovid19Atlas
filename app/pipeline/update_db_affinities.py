import pandas as pd

import sys
import os

script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append("/".join(script_dir.split("/")[:-2]))

from app import db
from app.api.models import Peptides, HLAAllelesPeptides

gisaid_id = sys.argv[1]

for hla_class in ["HLA-I"]:
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
