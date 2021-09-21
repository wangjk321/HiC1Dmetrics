juicer=$1 #juicertools path
hic=$2 #.hic file
chr=$3 #number (1,2,3....)
res=$4 #resolution
gt=$5 #genome_table
outname=$6

#step1 run hiccups
java -Xms512m -Xmx20480m -jar $juicer hiccups --cpu $hic looptemp -c $chr --threads 8

#step2 build fragment file for fit-hic
chromosome=chr$chr
length=`awk '$1 == "'$chromosome'" {print $2}' $gt`
endlength=`expr $length + $res`
mid=`expr $res / 2`

seq 0 $res $length > start
seq $mid $res $endlength > end
paste start end |awk -v ch=$chr '{printf "%s\t%d\t%d\t%s\t%s\n", ch,$1,$2,1,1}' |\
    gzip > fragment.temp.gz

seq $res $res $endlength > end2
paste start end2 |awk -v ch=$chromosome '{printf "%s\t%d\t%d\n", ch,$1,$2}'  > genome.split

#step3 make anchor.bed
touch Anchor.txt
cat looptemp/merged_loops.bedpe |sed '1d'|sed '1d'|cut -f 1-3 > Anchor.txt
cat looptemp/merged_loops.bedpe |sed '1d'|sed '1d'|cut -f 4-6 >> Anchor.txt
awk '{print "chr"$1"\t"$2"\t"$3}' Anchor.txt | sortBed > Anchor.bed

bedtools coverage -a genome.split -b Anchor.bed|cut -f 1-4 > $outname.bedGraph

rm -rf looptemp
rm end2 genome.split Anchor.bed Anchor.txt
rm start end fragment.temp.gz
