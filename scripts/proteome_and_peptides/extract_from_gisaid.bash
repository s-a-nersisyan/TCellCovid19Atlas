# Input:
# 1) GISAID Accesion ID (e.g., EPI_ISL_402125)
# 2) Output directory
# Output:
# 1) FASTA with protein sequences (proteome.fasta)
# 2) CSV table with HLA-I peptides (8-14 aminoacids, peptides_HLA-I.csv)
# 3) CSV table with HLA-II peptides (15-20 aminoacids, peptides_HLA-II.csv)

GISAID_ID=$1
OUT_DIR=$(realpath $2)

# Go to the script directory and load config file
cd $(dirname $0)
source ../config.bash

# Create output directory
mkdir -p $OUT_DIR

# Extract viral proteome from GISAID
if [ "$GISAID_ID" == "all" ]; then
    cat $GISAID_ALLPROT_PATH > $OUT_DIR/proteome.fasta
    # Generate 8-14-mers for HLA-I
    python3 generate_unique_peptides.py $OUT_DIR/proteome.fasta "8,9,10,11,12,13,14" > $OUT_DIR/peptides_HLA-I.csv
    # Generate 15-20-mers for HLA-II
    python3 generate_unique_peptides.py $OUT_DIR/proteome.fasta "15,16,17,18,19,20" > $OUT_DIR/peptides_HLA-II.csv

else
    grep -iF -A 1 "$GISAID_ID" $GISAID_ALLPROT_PATH > $OUT_DIR/proteome.fasta
    # Generate 8-14-mers for HLA-I
    python3 generate_peptides.py $OUT_DIR/proteome.fasta "8,9,10,11,12,13,14" > $OUT_DIR/peptides_HLA-I.csv
    # Generate 15-20-mers for HLA-II
    python3 generate_peptides.py $OUT_DIR/proteome.fasta "15,16,17,18,19,20" > $OUT_DIR/peptides_HLA-II.csv
fi
