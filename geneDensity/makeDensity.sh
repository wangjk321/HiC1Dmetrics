touch genometable.split
res=$1

for i in `seq 1 22` X
do
	chromosome=chr$i
	length=`awk '$1 == "'$chromosome'" {print $2}' genome_table`
	endlength=`expr $length + $res`
	
	seq 0 $res $length > start
	seq $res $res $endlength > end
	paste start end | sed "s/^/${chromosome}\t/g"  >> genometable.split
done

bedtools coverage -a genometable.split -b refFlat.bed |cut -f 1-4 > geneDensity$res.txt


rm start end genometable.split
