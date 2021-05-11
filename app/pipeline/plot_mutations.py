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
    ).filter(HLAAllelesPeptides.affinity <= AFFINITY_THRESHOLD)
    if (protein != None and protein != "Summary"):
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
        
        increase = increase / (tight_binders_num + 1) * 100
        decrease = decrease / (tight_binders_num + 1) * 100
        
        # if (tight_binders_num == 0):
        #     increase = decrease = 0
        # else:
        #     increase = increase / tight_binders_num * 100
        #     decrease = decrease / tight_binders_num * 100
        
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
 

    max_val.append(max(
        abs(np.max(df_prc["Weaker binding"])),
        abs(np.min(df_prc["Stronger binding"]))
    ))
   
    df_prc.plot.barh(
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
    
    # if (protein == None):
    #     protein = "Summary"
    
    plt.savefig(
        "../static/plots/{}/{}/{}.png".format(
            gisaid_id,
            hla_class,
            protein
        )
    )
    
    plt.close()


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
        HLAAllelesPeptides.hla_allele == allele
    ).filter(
        HLAAllelesPeptides.affinity <= AFFINITY_THRESHOLD   
    )
 
    df_num = pd.DataFrame(columns=["Protein", "Weaker binding", "Stronger binding"])
    df_prc = pd.DataFrame(columns=["Protein", "Weaker binding", "Stronger binding"])

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
            (report["Allele"] == protein) &
            (report["Ref aff"] > AFFINITY_THRESHOLD) &
            (report["Mut aff"] <= AFFINITY_THRESHOLD)
        ])
        
        # if (tight_binders_num == 0):
        #     increase = decrease = 0
        # else:
        #     increase = increase / tight_binders_num * 100
        #     decrease = decrease / tight_binders_num * 100
        
        increase = increase / (tight_binders_num + 1) * 100
        decrease = decrease / (tight_binders_num + 1) * 100
    

        df_prc.loc[len(df_prc)] = [protein, increase, decrease]
     
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
    
    df_num = df_num.set_index(df_num["Protein"])
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


    max_val.append(max(
        abs(np.max(df_prc["Weaker binding"])),
        abs(np.min(df_prc["Stronger binding"]))
    ))
 
    df_prc.plot.barh(
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

