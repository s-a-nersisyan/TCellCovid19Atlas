# Input:
# 1) Path to CSV table with HLA-I peptides
# 2) Number of processes
# 3) Output directory
# Output:
# CSV table with predicted affinities (binding_affinities_HLA-I.csv)

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
	"$netMHCpan -BA -p $OUT_DIR/peptides.txt -a {} > $OUT_DIR/raw_predictions/{}.txt"

python3 aggregate_HLA-I_affinities.py $OUT_DIR/raw_predictions/ > $OUT_DIR/binding_affinities_HLA-I.csv
rm -rf $OUT_DIR/raw_predictions
