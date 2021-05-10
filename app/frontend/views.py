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
    
    # Load report
    df = pd.read_csv("{}/{}/report.csv".format(app.config["PIPELINE_PATH"], gisaid_id))
    mut_df = pd.read_csv("{}/{}/mutations.csv".format(app.config["PIPELINE_PATH"], gisaid_id))
     
    alleles = sorted(set(df["Allele"]))
    if allele:
        df = df[df["Allele"] == allele]
    
    df["Ref peptide"] = df["Ref peptide"].fillna("-")
    df["Mut peptide"] = df["Mut peptide"].fillna("-")
    
    df["Ref aff"] = df["Ref aff"].fillna("-")
    df["Mut aff"] = df["Mut aff"].fillna("-")
    
    df["Ref aff"] = [int(aff) if type(aff) == float else aff for aff in df["Ref aff"]]
    df["Mut aff"] = [int(aff) if type(aff) == float else aff for aff in df["Mut aff"]]
    
    all_proteins = list(pd.read_csv("{}/{}/proteins.csv".format(app.config["PIPELINE_PATH"], gisaid_id))["Protein"])
    affected_proteins = set(df["Protein"])
    affected_proteins = [p for p in all_proteins if p in affected_proteins]
   
    return render_template(
        "variant_stats.html",
        df=df, mut_df=mut_df, gisaid_id=gisaid_id,
        all_proteins=all_proteins, affected_proteins=affected_proteins,
        alleles=alleles, allele=allele
    )

@frontend.route("/<gisaid_id>/download_report", methods=["GET"])
def download_report(gisaid_id):
    
    # Load report
    with open("{}/{}/report.csv".format(app.config["PIPELINE_PATH"], gisaid_id)) as rf:
        report = rf.read()
    
    return Response(
        report,
        mimetype="text/csv",
        headers={
            "Content-disposition":
            "attachment; filename={}.csv".format(gisaid_id)
        }
    )

