from . import frontend
from app import app

from flask import render_template
from flask import request
from flask import Response

import pandas as pd
import os
import json

@frontend.route("/", methods=["GET"])
def show_index():
    # Load table with available COVID-19 lineages
    lineages = pd.read_csv("{}/lineages.csv".format(app.config["PIPELINE_PATH"]))
    
    return render_template("index.html", lineages=lineages)


@frontend.route("/<gisaid_id>", methods=["GET"])
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

@frontend.route("/<gisaid_id>/download", methods=["GET"])
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
    
    
