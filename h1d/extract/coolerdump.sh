#!/bin/bash

cool=$1
binsize=$2
chr=$3
gt=$4
outname=$5

chrlen=$(cooler dump -t chroms $cool |head -1 |cut -f 1 |awk '{print length}')

pwd=$(cd $(dirname $0) && pwd)
chrnum=$(echo $chr |sed 's/chr//g')

if [ $chrlen > 3 ]
then
  cooler dump -r $chr --join $cool | cut -f 2,5,7  > sparse.temp
else
  cooler dump -r $chrnum --join $cool | cut -f 2,5,7  > sparse.temp
fi

mkdir -p $outname/${binsize}

outdir=$outname/${binsize}/${chr}.matrix.gz

$pwd/convert_JuicerDump_to_dense.py sparse.temp $outdir $gt $chr $binsize

rm sparse.temp
