from . import frontend
from app import app

from flask import render_template

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
    # Load report
    df = pd.read_csv("{}/{}/report.csv".format(app.config["PIPELINE_PATH"], gisaid_id))
    
    df["Ref peptide"] = df["Ref peptide"].fillna("-")
    df["Mut peptide"] = df["Mut peptide"].fillna("-")
    
    df["Ref aff"] = df["Ref aff"].fillna("-")
    df["Mut aff"] = df["Mut aff"].fillna("-")
    
    df["Ref aff"] = [int(aff) if type(aff) == float else aff for aff in df["Ref aff"]]
    df["Mut aff"] = [int(aff) if type(aff) == float else aff for aff in df["Mut aff"]]

    return render_template("variant_stats.html", df=df)
