import pandas as pd
import sys
import os

ref_diff_tb_path = sys.argv[1]
mut_diff_tb_path = sys.argv[2]
output_dir = sys.argv[3]

ref_diff_tb = pd.read_csv(ref_diff_tb_path, sep=",")
mut_diff_tb = pd.read_csv(mut_diff_tb_path, sep=",")

# swap mut and ref peptides and indexes for mut table
mut_diff_tb = mut_diff_tb.rename(columns={
    "ref_pep": "mut_pep",
    "mut_pep": "ref_pep",
    "ref_aln_pep": "mut_aln_pep",
    "mut_aln_pep": "ref_aln_pep",
    "ref_start": "mut_start",
    "ref_end": "mut_end",
    "mut_start": "ref_start",
    "mut_end": "ref_end",
})

output_tb = pd.concat(
    [ref_diff_tb, mut_diff_tb],
    ignore_index=True
)

output_tb = output_tb.drop_duplicates(
    subset=["ref_pep"]
)

output_tb.to_csv(os.path.join(output_dir, "diff.csv"),
                 sep=",", index=None)
