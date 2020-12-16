res=$1

touch genometable.split

for i in `seq 1 19` X
do
	chromosome=chr$i
	length=`awk '$1 == "'$chromosome'" {print $2}' mm10_genome_table`
	endlength=`expr $length + $res`
	
	seq 0 $res $length > start
	seq $res $res $endlength > end
	paste start end | sed "s/^/${chromosome}\t/g"  >> genometable.split
done

bedtools coverage -a genometable.split -b mm10_refFlat.bed |cut -f 1-4 > mm10_geneDensity$res.txt


rm start end genometable.split
