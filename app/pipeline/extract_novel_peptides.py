import pandas as pd

import sys
import os

script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append("/".join(script_dir.split("/")[:-2]))

from app import db
from app.api.models import Peptides

gisaid_id = sys.argv[1]
peptide_len_range = {
    "HLA-I": [8, 14],
    "HLA-II": [15, 20]
}

for hla_class in ["HLA-I", "HLA-II"]:
    df = pd.read_csv(os.path.join(script_dir, gisaid_id, "diff_{}.csv".format(hla_class)))

    novel_peptides = df["mut_pep"].unique()
    novel_peptides = {
        peptide for peptide in novel_peptides
        if not pd.isna(peptide) and peptide_len_range[hla_class][0] <= len(peptide) <= peptide_len_range[hla_class][1]
    }

    existing_peptides = {peptide.peptide for peptide in Peptides.query.all()}
    novel_peptides -= existing_peptides

    novel_peptides = ["Peptide"] + sorted(list(novel_peptides))
    f_out = open(
        "{}/novel_peptides_{}.csv".format(
            os.path.join(script_dir, gisaid_id), hla_class
        ), "w"
    )
    f_out.write("\n".join(novel_peptides))
    f_out.close()
