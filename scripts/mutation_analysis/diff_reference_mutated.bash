REF_DIR=$1
MUT_DIR=$2
OUT_DIR=$(realpath $2)

cd $(dirname $0)

HLA_classes="HLA-I HLA-II"

for HLA_class in $HLA_classes; do
	mkdir -p $OUT_DIR/ref_mut
	python3 get_diff_peptides.py \
		$REF_DIR/proteome.fasta $MUT_DIR/proteome.fasta \
		$REF_DIR/peptides_${HLA_class}.csv $MUT_DIR/peptides_${HLA_class}.csv \
		$OUT_DIR/ref_mut

	mkdir -p $OUT_DIR/mut_ref
	python3 get_diff_peptides.py \
		$MUT_DIR/proteome.fasta $REF_DIR/proteome.fasta \
		$MUT_DIR/peptides_${HLA_class}.csv $REF_DIR/peptides_${HLA_class}.csv \
		$OUT_DIR/mut_ref

	python3 merge_diff_tables.py $OUT_DIR/ref_mut/diff.csv $OUT_DIR/mut_ref/diff.csv $OUT_DIR
	mv $OUT_DIR/ref_mut/ref_aln.fasta $OUT_DIR/ref_aln.fasta
	mv $OUT_DIR/ref_mut/mut_aln.fasta $OUT_DIR/mut_aln.fasta
	mv $OUT_DIR/diff.csv $OUT_DIR/diff_${HLA_class}.csv

	rm -rf $OUT_DIR/ref_mut
	rm -rf $OUT_DIR/mut_ref
done
