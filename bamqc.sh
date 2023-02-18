#!/bin/bash

bedfile=$1
bamfile=$2
sample=$3

# mapq: 20; baseq: 20; rmdup: true; overlaps: true; indels: false; flank: 250

genome=/home/chenshulin/Project/Data/gatk-bundle-hg19/ucsc.hg19.fasta
covs="1,50,200,500"

mkdir bamqc
# picard
tmpname=`basename $bedfile`
tmpname=${tmpname%.*}
if [ ! -e bamqc/${sample}_hsMetrics.txt ]; then
  gatk BedToIntervalList -I $bedfile -O bamqc/${tmpname}.interval_list -SD ${genome%.*}.dict 
  gatk CollectHsMetrics -I $bamfile -R $genome -O bamqc/${sample}_hsMetrics.txt \
    --BAIT_INTERVALS bamqc/${tmpname}.interval_list \
    --TARGET_INTERVALS bamqc/${tmpname}.interval_list \
    --CLIP_OVERLAPPING_READS true \
    --INCLUDE_INDELS false \
    --MINIMUM_BASE_QUALITY 20 \
    --MINIMUM_MAPPING_QUALITY 20 \
    --NEAR_DISTANCE 250 \
    --PER_BASE_COVERAGE bamqc/${sample}_hsBase.txt \
    --PER_TARGET_COVERAGE bamqc/${sample}_hsTarget.txt
fi

if [ ! -e bamqc/${sample}_libComplex.txt ]; then
  gatk EstimateLibraryComplexity -I $bamfile -O bamqc/${sample}_libComplex.txt
fi

if [ ! -e bamqc/${sample}_insertSize.txt ]; then
  gatk CollectInsertSizeMetrics -I $bamfile -O bamqc/${sample}_insertSize.txt -H bamqc/${sample}_insertSize.pdf
fi

# bamdst
if [ ! -e bamqc/depth.tsv.gz ]; then
  bamdst -q 20 -f 250 -p $bedfile -o bamqc $bamfile
fi

# mosdpth
if [ ! -e bamqc/${sample}.thresholds.bed.gz ]; then
  mosdepth -b $bedfile -x -Q 20 -T ${covs} bamqc/${sample} $bamfile
fi


