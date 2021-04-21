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
from app.api.models import PeptidesPositions

import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams["font.sans-serif"] = "Arial"
mpl.rcParams["figure.dpi"] = 300

AFFINITY_THRESHOLD = 50

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

tight_binders = db.session.query(
    PeptidesPositions,
    HLAAllelesPeptides
).filter(HLAAllelesPeptides.affinity <= AFFINITY_THRESHOLD)
if (protein != None):
    tight_binders = tight_binders.filter(
        PeptidesPositions.protein == protein
    )


df_num = pd.DataFrame(columns=["Allele", "Weaker binding", "Stronger binding"])
df_prc = pd.DataFrame(columns=["Allele", "Weaker binding", "Stronger binding"])
for allele in sorted(alleles):
    tight_binders_num = tight_binders.filter(
        HLAAllelesPeptides.hla_allele == allele
    ).filter(PeptidesPositions.peptide == HLAAllelesPeptides.peptide).count()
    
    increase = len(report[
        (report["Allele"] == allele) &
        (report["Ref aff"] <= AFFINITY_THRESHOLD) &
        (report["Mut aff"] > AFFINITY_THRESHOLD)
    ])

    decrease = len(report[
        (report["Allele"] == allele) &
        (report["Ref aff"] > AFFINITY_THRESHOLD) &
        (report["Mut aff"] <= AFFINITY_THRESHOLD)
    ])
    
    if (tight_binders_num == 0):
        increase = decrease = 0
    else:
        increase = increase / tight_binders_num * 100
        decrease = decrease / tight_binders_num * 100

    df_prc.loc[len(df_prc)] = [allele, increase, decrease]


    increase = len(report[
        (report["Allele"] == allele) &
        (report["Ref aff"] < report["Mut aff"])
    ])

    decrease = len(report[
        (report["Allele"] == allele) &
        (report["Ref aff"] > report["Mut aff"])
    ])

    df_num.loc[len(df_num)] = [allele, increase, decrease]


df_prc = df_prc.set_index(df_prc["Allele"])
# df_prc = df_prc.sort_values(by=["Weaker binding", "Stronger binding"], ascending=True)
df_prc["Stronger binding"] = -df_prc["Stronger binding"]

df_num = df_num.set_index(df_num["Allele"])
df_num = df_num.sort_values(by=["Weaker binding", "Stronger binding"], ascending=True)
max_val = max(
    np.max(df_num["Weaker binding"]),
    np.max(df_num["Stronger binding"])
) + 1
df_num["Stronger binding"] = -df_num["Stronger binding"]
df_prc = df_prc.loc[df_num.index]


sns.set()
fig, axes = plt.subplots(nrows=1, ncols=2)
df_prc.plot.barh(
    stacked=True,
    ylabel="Percentage of tight-binding peptides", xlabel="",
    figsize=(12, len(df_prc) / 4),
    ax=axes[1],
    xlim=[-100, 100]
)
axes[1].legend(loc='lower right')
axes[1].set_xlabel("Percentage of tight-binding peptides")

df_num.plot.barh(
    stacked=True,
    ylabel="Number of peptides", xlabel="",
    figsize=(12, len(df_num) / 4),
    ax=axes[0],
    xlim=[-max_val, max_val],
    legend=False
)
# axes[0].legend(loc='lower right')
axes[0].set_xlabel("Number of peptides")

# if (protein != None):
#     plt.title(protein)
plt.tight_layout()

if (protein != None):
    protein += "_"
else:
    protein = ""

plt.savefig(
    os.path.join("{}/plots/{}report_summary.png".format(
        gisaid_id,
        protein
    ))
)
plt.close()

