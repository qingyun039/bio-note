#!/bin/bash

# -c 32 --mem 40G

seqdir=$1
clip=$2
lane=$3
barcode=$4
sample=$5

genome=/home/chenshulin/Project/Data/gatk-bundle-hg19/ucsc.hg19.fasta

if [ $lane == 'merge' ]
then
mkdir merge
cat ${seqdir}/${clip}/L0?/${clip}_L0?_${barcode}_1.fq.gz > merge/${clip}_${barcode}_1.fq.gz
cat ${seqdir}/${clip}/L0?/${clip}_L0?_${barcode}_2.fq.gz > merge/${clip}_${barcode}_2.fq.gz
read1=merge/${clip}_${barcode}_1.fq.gz
read2=merge/${clip}_${barcode}_2.fq.gz
else
read1=${seqdir}/${clip}/${lane}/${clip}_${lane}_${barcode}_1.fq.gz
read2=${seqdir}/${clip}/${lane}/${clip}_${lane}_${barcode}_2.fq.gz
fi

echo $read1 $read2

mkdir align
fastp \
  --thread 6 \
  -i ${read1} \
  -I ${read2} \
  --stdout \
  -j align/${sample}.json \
  -h align/${sample}.html \
  2> ${sample}.log \
| bwa-mem2 mem -M -Y -t 32 -p \
  -R "@RG\tID:${sample}\tPL:MGI\tLB:${clip}\tSM:${sample}" \
  $genome - 2>> ${sample}.log \
| samblaster 2>> ${sample}.log \
| samtools sort --write-index -@ 8 -O BAM -o align/${sample}.bam##idx##align/${sample}.bam.bai 2>> ${sample}.log
