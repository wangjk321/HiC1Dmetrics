#!/bin/bash

cool=$1
binsize=$2
chr=$3
gt=$4
outname=$5

pwd=$(cd $(dirname $0) && pwd)

echo Start dump matrix
#chrnum=$(echo $chr |sed 's/chr//g')
cooler dump -r $chr --join $cool | cut -f 2,5,7  > sparse.temp
#|| cooler dump -r $chrnum --join $cool | cut -f 2,5,7 > sparse.temp
echo Finish dump matrix

mkdir -p $outname/${binsize}

outdir=$outname/${binsize}/${chr}.matrix.gz

$pwd/convert_JuicerDump_to_dense.py sparse.temp $outdir $gt $chr $binsize

rm sparse.temp
