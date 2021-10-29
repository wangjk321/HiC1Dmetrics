#!/bin/bash

cool=$1
binsize=$2
chr=$3
gt=$4
outname=$5


pwd=$(cd $(dirname $0) && pwd)
chrnum=$(echo $chr |sed 's/chr//g')
chrlen=$(cooler dump -t chroms $cool |head -1 |cut -f 1 |awk '{print length}')
echo $chrlen

if [ $chrlen -gt 3 ]
then
  echo 'chang'
  cooler dump -r $chr --join $cool | cut -f 2,5,7  > sparse.temp
else
  echo 'duan'
  cooler dump -r $chrnum --join $cool | cut -f 2,5,7  > sparse.temp
fi

mkdir -p $outname/${binsize}

outdir=$outname/${binsize}/${chr}.matrix.gz

$pwd/convert_JuicerDump_to_dense.py sparse.temp $outdir $gt $chr $binsize

rm sparse.temp
