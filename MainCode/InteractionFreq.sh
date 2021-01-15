juicer=$1 #juicertools path
hic=$2 #.hic file
chr=$3 #number (1,2,3....)
res=$4 #resolution
gt=$5 #genome_table
outname=$6

#step1 dump interaction from .hic
java -Xms512m -Xmx20480m -jar $juicer dump observed KR $hic $chr $chr BP $res dump.temp.txt

#step2 convert hic dump file to fit-Hi-C input
cat dump.temp.txt | sed 's/NaN/0/g'| \
    awk -v chr1=$chr -v chr2=$chr '{ printf "%s\t%s\t%s\t%s\t%s\n", chr1,$1,chr2,$2,$3}' |\
    gzip > fithic.temp.gz

#step3 build fragment file for fit-hic
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

#step4 run Fit-Hi-C
fithic -i fithic.temp.gz -f fragment.temp.gz -r 50000 -o .

#step5 extract significant interactions
touch Anchor.txt
zcat < FitHiC.spline_pass1.res$res.significances.txt.gz|awk '$7<0.05{print $0}'|cut -f 2 > Anchor.txt
zcat < FitHiC.spline_pass1.res$res.significances.txt.gz|awk '$7<0.05{print $0}'|cut -f 4 >> Anchor.txt

awk '{print "'$chromosome'""\t"$1"\t"$1+"'$res'"}' Anchor.txt > Anchor.bed
bedtools coverage -a genome.split -b Anchor.bed|cut -f 1-4 > $outname.bedgraph

rm dump.temp.txt start end fragment.temp.gz fithic.temp.gz
rm FitHiC.fithic* FitHiC.spline_pass1.*.significances.txt.gz
rm end2 genome.split Anchor.bed Anchor.txt

#sh InteractionFreq.sh /Users/wangjiankang/Documents/localrun/juicer_tools_1.19.02.jar  /Users/wangjiankang/Documents/localrun/MCF7_Ctrl.hic 21 50000 /Users/wangjiankang/Documents/localrun/genome_table outname
