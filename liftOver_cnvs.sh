#!/bin/bash

python get_pdf_tables.py

for f in `ls *_array_cnvs.csv`; do
    sample=`echo $f | cut -d'_' -f1`
    echo starting $sample
    awk -F"," 'NR>1 {OFS="\t";print "chr"$2,$4,$5,$1}' $f > $sample.bed
    ./liftOver -multiple $sample.bed hg19ToHg38.over.chain.gz ${sample}_grch38.bed ${sample}_unmapped.bed
    exNum=`grep $sample sw_families.tsv | cut -f1`
done

python merge_builds.py
