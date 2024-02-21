#!/bin/bash
cmdname=`basename $0`
function usage()
{
    echo "$cmdname <norm> <matrixdir> <hic file> <binsize> <build>" 1>&2
}

if [ $# -ne 8 ]; then
  usage
  exit 1
fi

norm=$1
matrixdir=$2
hic=$3
binsize=$4
build=$5
juicer=$6
chrlist=$7
outname=$8
randomstr=`date +%s%N | md5sum | head -c 8 #`
uniqID=$chrlist$randomstr

pwd=$(cd $(dirname $0) && pwd)
juicertool="java -Xms512m -Xmx20480m -jar $juicer"
gt=$build
#chrlist=`cut -f 1 $gt`
#chrlist=$(/work/git/script_rnakato/getchr_from_genometable.sh $gt)
length=`awk '$1 == "'$chrlist'" {print $2}' $gt`
seq 0 $binsize $length > rowname.$uniqID.temp
cat rowname.$uniqID.temp|  tr '\n' '\t'|sed "s/^/\t/g"| sed "s/$/\n/g" > colname.$uniqID.temp
#sed '1{x;p;x}' rowname.temp > rowname2.temp

dir=$matrixdir/${outname}/$binsize
mkdir -p $dir
for chr in $chrlist #$(seq 1 $chrnum) X
do
    if test $chr = "chrY" -o $chr = "chrM" -o $chr = "chrMT" ;then continue; fi

    chr=$(echo $chr | sed -e 's/chr//g')
    echo "chr$chr"
    for type in observed #oe
    do
	     tempfile=$dir/$type.$norm.chr$chr.txt

       $juicertool dump $type $norm $hic $chr $chr BP $binsize $tempfile -d
       if [ $? -ne 0 ]; then
       #if ! test -s $tempfile; then
         echo try to use "chr1,chr2..." instead of "1,2..."
         # 如果执行错误，则执行第二个命令
         $juicertool dump $type $norm $hic chr$chr chr$chr BP $binsize $tempfile -d
       fi

       if [ $? -ne 0 ]; then
         echo "Are you using the correct version of Juicertools ??"
       fi

       if test -s $tempfile; then
		       paste rowname.$uniqID.temp $tempfile > pluscol.$uniqID.temp
		       cat colname.$uniqID.temp pluscol.$uniqID.temp |sed 's/\t$//g'| gzip > $dir/$type.$norm.chr$chr.matrix.gz
            #$pwd/convert_JuicerDump_to_dense.py \
	           #	$tempfile $dir/$type.$norm.chr$chr.matrix.gz $gt chr$chr $binsize
      fi
        rm rowname.$uniqID.temp pluscol.$uniqID.temp colname.$uniqID.temp $tempfile
    done
#    for type in expected norm
#    do
#        $juicertool dump $type $norm $hic.hic $chr BP $binsize $dir/$type.$norm.chr$chr.matrix -d
#    done
done
