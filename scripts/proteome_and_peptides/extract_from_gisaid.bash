# Input:
# 1) Path to GISAID "allprot" file
# 2) GISAID Accesion ID (e.g., EPI_ISL_402125)
# 3) Output directory
# Output:
# 1) FASTA with protein sequences (proteome.fasta)
# 2) CSV table with MHC-I peptides (8-14 aminoacids, peptides_MHC-I.csv)
# 3) CSV table with MHC-II peptides (15-20 aminoacids, peptides_MHC-II.csv)

GISAID_ALLPROT_PATH=$1
GISAID_ID=$2
OUT_DIR=$3

mkdir -p $OUT_DIR

# Extract viral proteome from GISAID
grep -iF -A 1 "$GISAID_ID" $GISAID_ALLPROT_PATH > $OUT_DIR/proteome.fasta

# Generate 8-14-mers for MHC-I
python3 generate_peptides.py $OUT_DIR/proteome.fasta "8,9,10,11,12,13,14" > $OUT_DIR/peptides_MHC-I.csv
# Generate 15-20-mers for MHC-II
python3 generate_peptides.py $OUT_DIR/proteome.fasta "15,16,17,18,19,20" > $OUT_DIR/peptides_MHC-II.csv
