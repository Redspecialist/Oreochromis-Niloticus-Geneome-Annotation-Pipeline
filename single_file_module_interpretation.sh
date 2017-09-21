INPUT_1="$1"
FILE_1=$(basename $INPUT_1)
FILE_1="${FILE_1%.*}"
DIR=$"./OUTPUT/"


ruby ./build_reference/convert_gene_names.rb $INPUT_1 "$DIR$FILE_1"
REF_1="$DIR$FILE_1.reference.txt"
echo $REF_1 "completed"

perl ./go_isolation/group_genes.pl $REF_1 "$DIR$FILE_1"
echo $REF_1 "go term tables completed"

/lustre/cichlid-labs/^C/R-3.1.0/bin/Rscript --no-restore-data count_find_module_genes.R $REF1 > "$DIR$FILE_1.module_count.txt"
echo $REF_1 "module count completed"
