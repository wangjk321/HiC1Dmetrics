res=$1
ref=$2
gt=$3
outname=$4

touch genometable.split

for i in `seq 1 22` X
do
	chromosome=chr$i
	length=`awk '$1 == "'$chromosome'" {print $2}' $gt`
	endlength=`expr $length + $res`
	
	seq 0 $res $length > start
	seq $res $res $endlength > end
	paste start end | sed "s/^/${chromosome}\t/g"  >> genometable.split
done

bedtools coverage -a genometable.split -b $ref |cut -f 1-4 > ${outname}.txt


rm start end genometable.split
