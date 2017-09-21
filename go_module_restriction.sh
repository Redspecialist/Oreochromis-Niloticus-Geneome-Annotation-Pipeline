OUTPUT_PREFIX="$3"
INPUT_1="$1"
INPUT_2=$2
FILE_1=$(basename $INPUT_1)
FILE_2=$(basename $INPUT_2)
FILE_1="${FILE_1%.*}"
FILE_2="${FILE_2%.*}"
DIR=$"./OUTPUT/"

REF_1="$DIR$FILE_1.reference.txt"
if [ ! -f $REF_1 ]
then
    ruby ./build_reference/convert_gene_names.rb $INPUT_1 "$DIR$FILE_1"
    echo $REF_1 "completed"
    ruby qtl_matching/match_qtls.rb qtl_matching/interestingQTLs.txt $REF_1 "$DIR$FILE_1.reference.QTL_matched.txt"
fi
REF_2="$DIR$FILE_2.reference.txt"

if [ ! -f $REF_2 ];
then
    ruby ./build_reference/convert_gene_names.rb $INPUT_2 "$DIR$FILE_2"
    echo $REF_2 "completed"
    ruby qtl_matching/match_qtls.rb qtl_matching/interestingQTLs.txt $REF_2 "$DIR$FILE_2.reference.QTL_matched.txt"
fi


perl ./module_comparison/compare_modules.pl $REF_1 $REF_2 $DIR$OUTPUT_PREFIX
echo $REF_1 $REF_2 "comparison completed"
ruby qtl_matching/match_qtls.rb qtl_matching/interestingQTLs.txt "$DIR$OUTPUT_PREFIX.comparison_sub_modules.out" "$DIR$OUTPUT_PREFIX.comparison_sub_modules.QTL_matched.out"


REF_1_filtered="$DIR$FILE_1.filtered.reference.txt"
if [ ! -f $REF_1_filtered ];
then
    perl ./go_isolation/group_genes.pl $REF_1 "$DIR$FILE_1"
    echo $REF_1 "complex tables completed"
    ruby qtl_matching/match_qtls.rb qtl_matching/interestingQTLs.txt $REF_1_filtered "$DIR$FILE_1.filtered.reference.QTL_matched.txt"
fi


REF_2_filtered="$DIR$FILE_2.filtered.reference.txt"
if [ ! -f $REF_2_filtered ];
then
    perl ./go_isolation/group_genes.pl $REF_2 "$DIR$FILE_2"
    echo $REF_2 "complex tables completed"
    ruby qtl_matching/match_qtls.rb qtl_matching/interestingQTLs.txt $REF_2_filtered "$DIR$FILE_2.filtered.reference.QTL_matched.txt"
fi

perl ./module_comparison/compare_modules.pl "$DIR$FILE_1.filtered.reference.txt" "$DIR$FILE_2.filtered.reference.txt" "$DIR$OUTPUT_PREFIX.filtered"
echo $REF_1 $REF_2 "filtered comparison completed"
ruby qtl_matching/match_qtls.rb qtl_matching/interestingQTLs.txt "$DIR$OUTPUT_PREFIX.filtered.comparison_sub_modules.out" "$DIR$OUTPUT_PREFIX.filtered.comparison_sub_modules.QTL_matched.out"
