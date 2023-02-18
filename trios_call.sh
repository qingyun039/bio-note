#!/bin/bash

bedfile=$1
fabam=$2
mobam=$3
chbam=$4

faid=$(basename $fabam)
faid=${faid%%.*}
moid=$(basename $mobam)
moid=${moid%%.*}
chid=$(basename $chbam)
chid=${chid%%.*}

echo -e "fam\t$faid\t0\t0\t1\t-9" > family.ped
echo -e "fam\t$moid\t0\t0\t2\t-9" >> family.ped
echo -e "fam\t$chid\t$faid\t$moid\t0\t-9" >> family.ped

genome=/home/chenshulin/Project/Data/gatk-bundle-hg19/ucsc.hg19.fasta
dbsnp=/home/chenshulin/Project/Data/gatk-bundle-hg19/dbsnp_138.hg19.vcf
gatk --java-options '-Xmx40g' HaplotypeCaller \
  -I ${fabam} \
  -I ${mobam} \
  -I ${chbam} \
  -O family.raw.vcf.gz \
  -L ${bedfile} \
  -R ${genome} \
  -D ${dbsnp} \
  -ped family.ped -ip 50

gatk --java-options '-Xmx40g' CalculateGenotypePosteriors \
  -V family.raw.vcf.gz \
  -O family.cgp.vcf.gz \
  -ped family.ped \
  --skip-population-priors
