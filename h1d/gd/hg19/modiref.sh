cat refFlat_full.txt |sed '1d'|awk -v OFS="\t" '{print $3,$5,$6,$2,$1,$4}' > refFlat_hg19.bed
