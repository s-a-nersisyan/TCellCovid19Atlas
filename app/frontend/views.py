from . import frontend
from app import app
from app.api.models import HLAAllelesPeptides
from sqlalchemy import or_

# from app.api.models import HLAAllelesPeptides

from flask import render_template
from flask import request
from flask import Response
from flask import jsonify

import pandas as pd
import pickle as pkl
from itertools import product
import os
import json


@frontend.route("/", methods=["GET"])
def show_index():
    # Load table with available COVID-19 lineages
    lineages = pd.read_csv("{}/lineages.csv".format(app.config["PIPELINE_PATH"]))
    
    return render_template("index.html", lineages=lineages)

@frontend.route("/stats/<gisaid_id>", methods=["GET"])
def show_report_page(gisaid_id):
    allele = request.args.get("allele")
    hla_class = request.args.get("hla_class")
    if not hla_class:
        hla_class = "HLA-I"

    # Load report
    df = pd.read_csv("{}/{}/report_{}.csv".format(app.config["PIPELINE_PATH"], gisaid_id, hla_class))
    mut_df = pd.read_csv("{}/{}/mutations.csv".format(app.config["PIPELINE_PATH"], gisaid_id))

    alleles = sorted(set(df["Allele"]))
    if allele:
        df = df[df["Allele"] == allele]
        allele_summary = pd.read_csv("{}/tables/{}/{}/Summary.csv".format(app.config["STATIC_PATH"], gisaid_id, hla_class))
        allele_summary = allele_summary[allele_summary["Allele"] == allele].to_dict(orient="records")[0]
    else:
        allele_summary=None

    df["Ref peptide"] = df["Ref peptide"].fillna("-")
    df["Mut peptide"] = df["Mut peptide"].fillna("-")
    
    df["Ref aff"] = df["Ref aff"].fillna("-")
    df["Mut aff"] = df["Mut aff"].fillna("-")
    
    df["Ref aff"] = [int(aff) if type(aff) == float else aff for aff in df["Ref aff"]]
    df["Mut aff"] = [int(aff) if type(aff) == float else aff for aff in df["Mut aff"]]
   
    all_proteins = list(pd.read_csv("{}/{}/proteins.csv".format(app.config["PIPELINE_PATH"], gisaid_id))["Protein"])
    affected_proteins = set(df["Protein"])
    affected_proteins = [p for p in all_proteins if p in affected_proteins]
    
    # Get variant name by gisaid id
    lineages = pd.read_csv("{}/lineages.csv".format(app.config["PIPELINE_PATH"]))
    lineage = lineages.loc[lineages["GISAID Accession ID"] == gisaid_id, "Lineage"].iloc[0]
   
    return render_template(
        "variant_stats.html",
        df=df, mut_df=mut_df, allele_summary=allele_summary, gisaid_id=gisaid_id, lineage=lineage,
        all_proteins=all_proteins, affected_proteins=affected_proteins,
        alleles=alleles, allele=allele, hla_class=hla_class
    )

@frontend.route("/stats/<gisaid_id>/download", methods=["GET"])
def download(gisaid_id):
    hla_class = request.args.get("hla_class")
    ident = request.args.get("type")
    
    # Load report
    with open("{}/{}/report_{}.csv".format(app.config["PIPELINE_PATH"], gisaid_id, hla_class)) as rf:
        report = rf.read()
    
    return Response(
        report,
        mimetype="text/csv",
        headers={
            "Content-disposition":
            "attachment; filename={}.csv".format(gisaid_id)
        }
    )
    
PROTEINS = [
    "All", "Spike", "N", "M", "E", "NS3", "NS6",
    "NS7a", "NS7b", "NS8", "NS9b", "NS9c",
    "NSP1", "NSP2", "NSP3", "NSP4", "NSP5",
    "NSP6", "NSP7", "NSP8", "NSP9", "NSP10",
    "NSP11", "NSP12", "NSP13", "NSP14", "NSP15", "NSP16"
]

def compare_variants(
    first_gisaid_id, second_gisaid_id,
    first_binders_dict, second_binders_dict,
    hla_alleles, protein
):
    gisaid_id_binders = {}
    iterator = zip(
        [first_gisaid_id, second_gisaid_id],
        [first_binders_dict, second_binders_dict]
    )

    for gisaid_id, binders_dict in iterator:
        _dict = binders_dict[protein] 

        binders = set()
        for allele in hla_alleles:
            if not (allele in _dict):
                continue
            
            binders.update(_dict[allele])
        
        gisaid_id_binders[gisaid_id] = binders
    
    first_score = len(gisaid_id_binders[first_gisaid_id])
    second_score = len(gisaid_id_binders[second_gisaid_id])
    
    return first_score, second_score

def hla_summary(
    first_gisaid_id, second_gisaid_id,
    hla_alleles, proteins
):
    summary = {}
    for protein in proteins:
        first_binders_dict = pkl.load(open("{}/{}/tight_binders.pkl".format(
            app.config["PIPELINE_PATH"], first_gisaid_id
        ), "rb"))
        
        second_binders_dict = pkl.load(open("{}/{}/tight_binders.pkl".format(
            app.config["PIPELINE_PATH"], second_gisaid_id
        ), "rb"))

        binder_analysis = {}
        for hla_allele in hla_alleles:
            binder_analysis[hla_allele] = compare_variants(
                first_gisaid_id, second_gisaid_id,
                first_binders_dict, second_binders_dict,
                [hla_allele], protein
            )

        binder_analysis["Summary"] = compare_variants(
            first_gisaid_id, second_gisaid_id,
            first_binders_dict, second_binders_dict,
            hla_alleles, protein
        )

        df = pd.DataFrame([
            {
                "Allele": hla_allele,
                "First": binder_analysis[hla_allele][0],
                "Second": binder_analysis[hla_allele][1]
            } for hla_allele in hla_alleles + ["Summary"]
        ])
        
        summary[protein] = df

    return summary

HLA_I_order = ["HLA-A", "HLA-B", "HLA-C"]
HLA_II_order = ["HLA-DRB1", "HLA-DQB1", "HLA-DPB1"]
def analyze_haplotype(alleles, hla_class):
    if hla_class == "HLA-I":
        hla_order = HLA_I_order 
    else:
        hla_order = HLA_II_order 
        
    hla_dict = {tp: [] for tp in hla_order}
    for allele in alleles:
        for tp in hla_order:
            if (
                (hla_class == "HLA-I" and allele.startswith(tp)) or
                (hla_class == "HLA-II" and tp.split("-")[-1] in allele)
            ):
                hla_dict[tp].append(allele)
                continue
    
    return hla_dict

def analyze_haplotype_frequency(alleles, hla_class):
    hla_dict = analyze_haplotype(alleles, hla_class)

    if hla_class == "HLA-I":
        hla_order = HLA_I_order 
    else:
        hla_order = HLA_II_order 
    
    for tp in hla_order:
        for i, allele in enumerate(hla_dict[tp]):
            hla_dict[tp][i] = allele.split("/")[-1].split("*")[-1]

        if not hla_dict[tp]:
            hla_dict[tp].append("")

        hla_dict[tp] = list(set(hla_dict[tp]))
    
    if hla_class == "HLA-I":
        frequency_df = pd.read_csv(
            "{}alleles/ABC_haplotype_region_frequency.tsv".format(app.config["STATIC_PATH"]),
            sep="\t"
        )
    else:
        frequency_df = pd.read_csv(
            "{}alleles/RQP_haplotype_region_frequency.tsv".format(app.config["STATIC_PATH"]),
            sep="\t"
        )
    
    for i, tp in enumerate(hla_order):
        frequency_df[tp] = frequency_df["Haplotype"].str.split("/").str[i]
    frequency_df = frequency_df.drop(columns=["Haplotype"])
    
    group_columns = [tp for tp in hla_order if hla_dict[tp][0]]

    haplotype_freq_df = pd.DataFrame(columns=frequency_df.columns) 
    
    hla_haplotypes = [hla_dict[tp] for tp in hla_order]
    for haplotype in product(*hla_haplotypes):
        to_append = frequency_df.loc[
            frequency_df[hla_order[0]].str.startswith(haplotype[0]) & \
            frequency_df[hla_order[1]].str.startswith(haplotype[1]) & \
            frequency_df[hla_order[2]].str.startswith(haplotype[2])
        ].groupby(group_columns).sum()
        
        for i, tp in enumerate(hla_order):
            to_append[tp] = haplotype[i]

        haplotype_freq_df = haplotype_freq_df.append(to_append, ignore_index=True)

    haplotype = "/".join([tp.split("-")[-1] for tp in hla_order])
    haplotype_freq_df[haplotype] = ""
    for i, tp in enumerate(hla_order):
        if len(haplotype_freq_df[tp]):
            if not haplotype_freq_df[tp].iloc[0]:
                haplotype_freq_df[haplotype] += "-"
            else:
                haplotype_freq_df[haplotype] += haplotype_freq_df[tp]

        if i < len(hla_order) - 1:
            haplotype_freq_df[haplotype] += "/"
    
    haplotype_freq_df = haplotype_freq_df.set_index(haplotype)
    haplotype_freq_df = haplotype_freq_df.drop(columns=hla_order)
    
    return haplotype_freq_df 
    
@frontend.route("/comparison/<first_gisaid_id>/<second_gisaid_id>", methods=["GET"])
def show_comparison_page(first_gisaid_id, second_gisaid_id):
    hla_alleles = request.args.get("hla_alleles") 
    hla_alleles = hla_alleles.split(",")
    hla_class = request.args.get("hla_class")
    
    known_alleles = set(json.load(open(
        "{}/alleles/alleles.json".format(app.config["STATIC_PATH"])
    )))

    hla_I_alleles = [
        hla for hla in hla_alleles if (not "D" in hla) and (hla in known_alleles)
    ]
    hla_II_alleles = [
        hla for hla in hla_alleles if ("D" in hla) and (hla in known_alleles)
    ]
 
    other_class = True

    if (not hla_I_alleles) and (not hla_II_alleles):
        other_class = False
        return "Unknown sequence of hla alleles", 400
    elif not hla_I_alleles:
        hla_class = "HLA-II"
        hla_alleles = hla_II_alleles
        other_class = False
    elif not hla_II_alleles:
        hla_class = "HLA-I"
        hla_alleles = hla_I_alleles
        other_class = False
    elif not hla_class:
        hla_class = "HLA-I"
        hla_alleles = hla_I_alleles
    elif hla_class == "HLA-II":
        hla_alleles = hla_II_alleles
    else:
        hla_class = "HLA-I"
        hla_alleles = hla_I_alleles

    analyze_haplotype_frequency(hla_alleles, hla_class)

    summary = hla_summary(
        first_gisaid_id, second_gisaid_id,
        hla_alleles, PROTEINS
    )

    lineages = pd.read_csv(
        "{}/lineages.csv".format(app.config["PIPELINE_PATH"])
    )

    haplotype_frequency = analyze_haplotype_frequency(
        hla_alleles, hla_class
    )
    
    return render_template(
        "variant_comparison.html",
        hla_summary=summary,
        proteins=PROTEINS,
        haplotype_frequency=haplotype_frequency,
        hla_alleles=request.args.get("hla_alleles"),
        hla_class=hla_class,
        first_gisaid_id=first_gisaid_id,
        second_gisaid_id=second_gisaid_id,
        lineages=lineages,
        other_class=other_class
    )

@frontend.route("/comparison/<first_gisaid_id>/<second_gisaid_id>/download", methods=["GET"])
def download_comparison(first_gisaid_id, second_gisaid_id):
    hla_alleles = request.args.get("hla_alleles") 
    hla_alleles = hla_alleles.split(",")
    protein = request.args.get("protein")

    if not protein in PROTEINS:
        protein = "All"

    gisaid_id_binders = {}
    df_peptides = []
    df_alleles = []
    for gisaid_id in [first_gisaid_id, second_gisaid_id]:
        _dict = pkl.load(open("{}/{}/tight_binders.pkl".format(
            app.config["PIPELINE_PATH"], gisaid_id
        ), "rb"))

        binders = set()
        for allele in hla_alleles:
            binders.update(_dict[protein][allele])
            
            df_peptides.extend(_dict[protein][allele])
            df_alleles.extend([allele] * len(_dict[protein][allele]))
        
        gisaid_id_binders[gisaid_id] = binders
    
    common_binders = gisaid_id_binders[first_gisaid_id] & \
        gisaid_id_binders[second_gisaid_id]
    
    disappeared_binders = gisaid_id_binders[first_gisaid_id] - \
        gisaid_id_binders[second_gisaid_id]
    
    appeared_binders = gisaid_id_binders[second_gisaid_id] - \
        gisaid_id_binders[first_gisaid_id]
    
    download_df = pd.DataFrame()
    download_df["Allele"] = df_alleles
    download_df["Peptide"] = df_peptides
    download_df = download_df.drop_duplicates(subset=["Allele", "Peptide"])

    df_type = []
    for _, row in download_df.iterrows():
        if row["Peptide"] in disappeared_binders:
            df_type.append("disappeared")
        elif row["Peptide"] in appeared_binders:
            df_type.append("appeared")
        else:
            df_type.append("common")
    
    download_df["Type"] = df_type

    db_filter = [
        (HLAAllelesPeptides.peptide==row["Peptide"]) &
        (HLAAllelesPeptides.hla_allele==row["Allele"]) \
        for _, row in download_df.iterrows()
    ]
    
    if len(download_df) > 0:
        download_df["Affinity"] = [
            row.affinity for row in HLAAllelesPeptides.query.filter(
                or_(*db_filter)
            ).all()
        ]
    else:
        download_df["Affinity"] = []

    return Response(
        download_df.to_csv(index=None, sep=","),
        mimetype="text/csv",
        headers={
            "Content-disposition":
            "attachment; filename={}_vs_{}_{}.csv".format(
                first_gisaid_id, second_gisaid_id,
                protein
            )
        }
    )

@frontend.route("/haplotypes", methods=["GET"])
def show_input_page():
    lineages = pd.read_csv(
        "{}/lineages.csv".format(app.config["PIPELINE_PATH"])
    )
    alleles = json.load(open(
        "{}/alleles/alleles.json".format(app.config["STATIC_PATH"])
    ))

    return render_template(
        "input_haplotypes.html",
        lineages=lineages,
        alleles=alleles
    )
