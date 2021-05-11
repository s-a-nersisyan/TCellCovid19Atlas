GISAID_ID=$1

# Go to the script directory
cd $(dirname $0)
OUT_DIR=$(pwd)/$GISAID_ID

# Extract peptides
../../scripts/proteome_and_peptides/extract_from_gisaid.bash $GISAID_ID $OUT_DIR

# Intersect peptides with reference ones
../../scripts/mutation_analysis/diff_reference_mutated.bash \
	../../data/reference_EPI_ISL_402125/ $OUT_DIR $OUT_DIR

# Extract novel peptides and predict HLA-I, HLA-II affinities
python3 extract_novel_peptides.py $GISAID_ID
../../scripts/HLA_peptide_binding/predict_HLA-I_affinities.bash \
	$GISAID_ID/novel_peptides_HLA-I.csv 1 $GISAID_ID

../../scripts/HLA_peptide_binding/predict_HLA-II_affinities.bash \
	$GISAID_ID/novel_peptides_HLA-II.csv 1 $GISAID_ID

# Write the results in DB
python3 update_db_affinities.py $GISAID_ID
python3 generate_report.py $GISAID_ID
