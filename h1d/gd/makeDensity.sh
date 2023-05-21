function usage(){
        echo "makeDensity.sh [-r resolution] [-g refSeq_gene] [-t genometable] [-o output_name]"
}

while getopts "r:g:t:o:" opt
do
        case $opt in
                r)
                        res=${OPTARG}
                        ;;
                g)
                        ref=${OPTARG}
                        ;;
		t)
			gt=${OPTARG}
			;;
		o)
			outname=${OPTARG}
			;;
                *)
                        usage
                        exit 1
                        ;;
        esac
done

if test $# -ne 8; then
    usage
    exit 0
fi

touch genometable.split
sed '1d' $ref| awk -v OFS="\t" '{print $3,$5,$6,$2,$1,$4}' >ref.bed

TAB="$(printf '\\\011')"

for i in `cat $gt | cut -f 1`
do
	chromosome=$i
	length=`awk '$1 == "'$chromosome'" {print $2}' $gt`
	endlength=`expr $length + $res`

	seq 0 $res $length > start
	seq $res $res $endlength > end
	paste start end | sed "s/^/${chromosome}${TAB}/g"  >> genometable.split
done

bedtools coverage -a genometable.split -b ref.bed |cut -f 1-4 > ${outname}.txt


rm start end ref.bed genometable.split
