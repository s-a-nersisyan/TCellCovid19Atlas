# Input:
# 1) Path to CSV table with MHC-I peptides
# 2) Number of processes
# 2) Output directory
# Output:
# 1) FASTA with protein sequences (proteome.fasta)
# 2) CSV table with MHC-I peptides (8-14 aminoacids, peptides_MHC-I.csv)
# 3) CSV table with MHC-II peptides (15-20 aminoacids, peptides_MHC-II.csv)

PEPTIDES_TABLE=$(realpath $1)
NUM_PROC=$2
OUT_DIR=$(realpath $3)

# Go to the script directory and load config file
cd $(dirname $0)
source ../config.bash

# Create output directory
mkdir -p $OUT_DIR
mkdir -p $OUT_DIR/raw_predictions

# Select peptides from the first column and remove header
cat $PEPTIDES_TABLE | cut -d "," -f 1 | tail -n +2 > $OUT_DIR/peptides.txt

cat $HLA_I_ALLELES_FILE | \
	parallel -j $NUM_PROC \
	"netMHCpan -BA -p $OUT_DIR/peptides.txt -a {} > $OUT_DIR/raw_predictions/{}.txt"

# TODO: merge results into a single table
