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

pwd=$(cd $(dirname $0) && pwd)
juicertool="java -Xms512m -Xmx20480m -jar $juicer"
gt=$build
#chrlist=`cut -f 1 $gt`
#chrlist=$(/work/git/script_rnakato/getchr_from_genometable.sh $gt)
length=`awk '$1 == "'$chrlist'" {print $2}' $gt`
seq 0 $binsize $length > rowname.temp
cat rowname.temp|  tr '\n' '\t'|sed "s/^/\t/g"| sed "s/$/\n/g" > colname.temp
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
	if test -s $tempfile; then
		paste rowname.temp $tempfile > pluscol.temp
		cat colname.temp pluscol.temp |sed 's/\t$//g'| gzip > $dir/$type.$norm.chr$chr.matrix.gz
            #$pwd/convert_JuicerDump_to_dense.py \
	#	$tempfile $dir/$type.$norm.chr$chr.matrix.gz $gt chr$chr $binsize
	fi
        rm *temp $tempfile
    done
#    for type in expected norm
#    do
#        $juicertool dump $type $norm $hic.hic $chr BP $binsize $dir/$type.$norm.chr$chr.matrix -d
#    done
done
