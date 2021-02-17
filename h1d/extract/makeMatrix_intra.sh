#!/bin/bash
cmdname=`basename $0`
function usage()
{
    echo "$cmdname <norm> <matrixdir> <hic file> <binsize> <build>" 1>&2
}

if [ $# -ne 5 ]; then
  usage
  exit 1
fi

norm=$1
matrixdir=$2
hic=$3
binsize=$4
build=$5

pwd=$(cd $(dirname $0) && pwd)
juicertool="java -Xms512m -Xmx20480m -jar /home/wang/software/juicer/juicer_tools_1.11.04_jcuda.0.8.jar"
gt=/work/Database/UCSC/$build/genome_table
chrlist=$(/work/git/script_rnakato/getchr_from_genometable.sh $gt)

dir=$matrixdir/Matrix/intrachromosomal/$binsize
mkdir -p $dir
for chr in $chrlist #$(seq 1 $chrnum) X
do
    if test $chr = "chrY" -o $chr = "chrM" -o $chr = "chrMT" ;then continue; fi

    chr=$(echo $chr | sed -e 's/chr//g')
    echo "chr$chr"
    for type in observed #oe
    do
	tempfile=$dir/$type.$norm.chr$chr.txt
        $juicertool dump $type $norm $hic $chr $chr BP $binsize $tempfile
	if test -s $tempfile; then
            $pwd/convert_JuicerDump_to_dense.py \
		$tempfile $dir/$type.$norm.chr$chr.matrix.gz $gt chr$chr $binsize
	fi
        rm $tempfile
    done
#    for type in expected norm
#    do
#        $juicertool dump $type $norm $hic.hic $chr BP $binsize $dir/$type.$norm.chr$chr.matrix -d
#    done
done
