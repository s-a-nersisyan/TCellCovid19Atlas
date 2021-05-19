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
REFERENCE_GISAID_ID = "EPI_ISL_402125"

# copied from generate_plot.py
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
# end of the block


def protein_plot(gisaid_id, hla_class, protein,):
    report = pd.read_csv("{}/report_{}.csv".format(
        gisaid_id,
        hla_class
    ))
    
    if (protein != "Summary"):
        report = report[report["Protein"] == protein]
    
    report = report.fillna(np.inf)
    
    alleles = set(report["Allele"])
    
    tight_binders = db.session.query(
        PeptidesPositions,
        HLAAllelesPeptides
    ).filter(
        PeptidesPositions.gisaid_id == REFERENCE_GISAID_ID
    ).filter(HLAAllelesPeptides.affinity <= AFFINITY_THRESHOLD)
    if (protein != None and protein != "Summary"):
        tight_binders = tight_binders.filter(
            PeptidesPositions.protein == protein
        )
    
    
    df_num = pd.DataFrame(columns=["Allele", "Weaker binding", "Stronger binding"])
    df_prc = pd.DataFrame(columns=["Allele", "Weaker binding", "Stronger binding", "Tight binding"])
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
        
        # increase = increase / (tight_binders_num + 1) * 100
        # decrease = decrease / (tight_binders_num + 1) * 100
        
        df_prc.loc[len(df_prc)] = [allele, increase, decrease, tight_binders_num] 
    
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
    df_prc["Stronger binding"] = -df_prc["Stronger binding"]
    
    df_num = df_num.set_index(df_num["Allele"])
    df_num = df_num.sort_values(by=["Weaker binding", "Stronger binding"], ascending=True)
    df_num["Stronger binding"] = -df_num["Stronger binding"]
    df_prc = df_prc.loc[df_num.index]
    
    sns.set()
    fig, axes = plt.subplots(nrows=1, ncols=2)
    max_val = []

    max_val.append(max(
        abs(np.max(df_num["Weaker binding"])),
        abs(np.min(df_num["Stronger binding"]))
    ))
   
    df_num.plot.barh(
        stacked=True,
        ylabel="Number of peptides", xlabel="",
        figsize=(12, 2 + len(df_num) / 4),
        ax=axes[0],
        legend=False
    )
    
    axes[0].set_xlabel("Number of peptides") 


    df_prc_draw = df_prc.copy()
    df_prc_draw["Weaker binding"] = df_prc_draw["Weaker binding"] /\
            (df_prc_draw["Tight binding"] + 1) * 100
    df_prc_draw["Stronger binding"] = df_prc_draw["Stronger binding"] /\
            (df_prc_draw["Tight binding"] + 1) * 100
    df_prc_draw = df_prc_draw[["Weaker binding", "Stronger binding"]]
    
    max_val.append(max(
        abs(np.max(df_prc_draw["Weaker binding"])),
        abs(np.min(df_prc_draw["Stronger binding"]))
    ))
 
    df_prc_draw.plot.barh(
        stacked=True,
        ylabel="Percentage of tight-binding peptides", xlabel="",
        figsize=(12, 2 + len(df_prc) / 4),
        ax=axes[1],
    )
    axes[1].legend(loc='lower right')
    axes[1].set_xlabel("Percentage of tight-binding peptides")
    
     
    # set ticks
    for j in [0, 1]:
        if (max_val[j] > 4):
            border = int(max_val[j] - 1) // 4 * 4 + 4
            step = border // 4
            axes[j].set_xticks(
                [step * i for i in range(-4, 5)]
            )
            axes[j].set_xlim(-border * 11 / 10, border * 11 / 10)
        else:
            axes[j].set_xticks(list(range(-int(max_val[j]), int(max_val[j]) + 1)))
            axes[j].set_xlim(-int(max_val[j]) - 1, int(max_val[j]) + 1)


    # set ticks' labels
    fig.canvas.draw()
    for i in range(len(axes)):
        labels = []
        for item in axes[i].get_xticklabels():
            if (item.get_text()[0] == "−"):
                labels.append(item.get_text()[1:])
            else:
                labels.append(item.get_text())
        axes[i].set_xticklabels(labels)
    
    plt.tight_layout()
    
    plt.savefig(
        "../static/plots/{}/{}/{}.png".format(
            gisaid_id,
            hla_class,
            protein
        )
    )
    
    plt.close()

    # save tables
    df_num["Stronger binding"] = -df_num["Stronger binding"]
    df_num = df_num.rename({
        "Weaker binding" : "Number of weaker binding",
        "Stronger binding" : "Number of stronger binding"
    }, axis="columns")
    
    df_prc["Stronger binding"] = -df_prc["Stronger binding"]
    df_prc = df_prc.rename({
        "Weaker binding" : "Number of weaker tight binding",
        "Stronger binding" : "Number of stronger tight binding",
        "Tight binding" : "Number of tight binding"
    }, axis="columns")
    
    df_num = df_num.drop(["Allele"], axis="columns")
    df_prc = df_prc.drop(["Allele"], axis="columns")
    df = df_num.merge(df_prc, left_index=True, right_index=True)
 
    df.to_csv("../static/tables/{}/{}/{}.csv".format(
        gisaid_id, hla_class, protein
    ))



def allele_plot(gisaid_id, hla_class, allele):
    report = pd.read_csv("{}/report_{}.csv".format(
        gisaid_id,
        hla_class
    ))
    
    report = report[report["Allele"] == allele]
    
    report = report.fillna(np.inf)
    
    proteins = set(report["Protein"])
    
    tight_binders = db.session.query(
        PeptidesPositions,
        HLAAllelesPeptides
    ).filter(
        PeptidesPositions.gisaid_id == REFERENCE_GISAID_ID
    ).filter(
        HLAAllelesPeptides.hla_allele == allele
    ).filter(
        HLAAllelesPeptides.affinity <= AFFINITY_THRESHOLD   
    )
 
    df_num = pd.DataFrame(columns=["Protein", "Weaker binding", "Stronger binding"])
    df_prc = pd.DataFrame(columns=["Protein", "Weaker binding", "Stronger binding", "Tight binding"])

    for protein in sorted(proteins):
        tight_binders_num = tight_binders.filter(
            PeptidesPositions.protein == protein
        ).filter(PeptidesPositions.peptide == HLAAllelesPeptides.peptide).count()
        
        increase = len(report[
            (report["Protein"] == protein) &
            (report["Ref aff"] <= AFFINITY_THRESHOLD) &
            (report["Mut aff"] > AFFINITY_THRESHOLD)
        ])
    
        decrease = len(report[
            (report["Protein"] == protein) &
            (report["Ref aff"] > AFFINITY_THRESHOLD) &
            (report["Mut aff"] <= AFFINITY_THRESHOLD)
        ])
        
        # increase = increase / (tight_binders_num + 1) * 100
        # decrease = decrease / (tight_binders_num + 1) * 100
    

        df_prc.loc[len(df_prc)] = [protein, increase, decrease, tight_binders_num]
     
        increase = len(report[
            (report["Protein"] == protein) &
            (report["Ref aff"] < report["Mut aff"])
        ])
    
        decrease = len(report[
            (report["Protein"] == protein) &
            (report["Ref aff"] > report["Mut aff"])
        ])
    
        df_num.loc[len(df_num)] = [protein, increase, decrease]
    
    
    df_prc = df_prc.set_index(df_prc["Protein"])
    df_prc["Stronger binding"] = -df_prc["Stronger binding"]
    
    df_num["Sort"] = [PROTEIN_SIGNIFICANCE[row["Protein"]] for _, row in df_num.iterrows()]
    df_num = df_num.sort_values("Sort", ascending="True")
    df_num = df_num.drop(columns=["Sort"])
    df_num = df_num.iloc[::-1] 

    df_num = df_num.set_index(df_num["Protein"])
    df_num["Stronger binding"] = -df_num["Stronger binding"]
    df_prc = df_prc.loc[df_num.index]
   
    
    sns.set()
    fig, axes = plt.subplots(nrows=1, ncols=2)
    max_val = []

    max_val.append(max(
        abs(np.max(df_num["Weaker binding"])),
        abs(np.min(df_num["Stronger binding"]))
    ))

    df_num.plot.barh(
        stacked=True,
        ylabel="Number of peptides", xlabel="",
        figsize=(12, 2 + len(df_num) / 4),
        ax=axes[0],
        legend=False
    )
    
    axes[0].set_xlabel("Number of peptides") 

        
    df_prc_draw = df_prc.copy()
    df_prc_draw["Weaker binding"] = df_prc_draw["Weaker binding"] /\
            (df_prc_draw["Tight binding"] + 1) * 100
    df_prc_draw["Stronger binding"] = df_prc_draw["Stronger binding"] /\
            (df_prc_draw["Tight binding"] + 1) * 100
    df_prc_draw = df_prc_draw[["Weaker binding", "Stronger binding"]]
    
    max_val.append(max(
        abs(np.max(df_prc_draw["Weaker binding"])),
        abs(np.min(df_prc_draw["Stronger binding"]))
    ))

    df_prc_draw.plot.barh(
        stacked=True,
        ylabel="Percentage of tight-binding peptides", xlabel="",
        figsize=(12, 2 + len(df_prc) / 4),
        ax=axes[1],
    )
    axes[1].legend(loc='lower right')
    axes[1].set_xlabel("Percentage of tight-binding peptides")

    # set ticks
    for j in [0, 1]:
        if (max_val[j] > 4):
            border = int(max_val[j] - 1) // 4 * 4 + 4
            step = border // 4
            axes[j].set_xticks(
                [step * i for i in range(-4, 5)]
            )
            axes[j].set_xlim(-border * 11 / 10, border * 11 / 10)
        else:
            axes[j].set_xticks(list(range(-int(max_val[j]), int(max_val[j]) + 1)))
            axes[j].set_xlim(-int(max_val[j]) - 1, int(max_val[j]) + 1)

    # set ticks' labels
    fig.canvas.draw()
    for i in range(len(axes)):
        labels = []
        for item in axes[i].get_xticklabels():
            if (item.get_text()[0] == "−"):
                labels.append(item.get_text()[1:])
            else:
                labels.append(item.get_text())
        axes[i].set_xticklabels(labels)
    
    plt.tight_layout()
        
    plt.savefig(
        "../static/plots/{}/{}/{}.png".format(
            gisaid_id,
            hla_class,
            allele
        )
    )
    
    plt.close()

    # save tables
    df_num["Stronger binding"] = -df_num["Stronger binding"]
    df_num = df_num.rename({
        "Weaker binding" : "Number of weaker binding",
        "Stronger binding" : "Number of stronger binding"
    }, axis="columns")
    
    df_prc["Stronger binding"] = -df_prc["Stronger binding"]
    df_prc = df_prc.rename({
        "Weaker binding" : "Number of weaker tight binding",
        "Stronger binding" : "Number of stronger tight binding",
        "Tight binding" : "Number of tight binding"
    }, axis="columns")
    
    df_num = df_num.drop(["Protein"], axis="columns")
    df_prc = df_prc.drop(["Protein"], axis="columns")
    df = df_num.merge(df_prc, left_index=True, right_index=True)
    
    df.to_csv("../static/tables/{}/{}/{}.csv".format(
        gisaid_id, hla_class, allele
    ))
