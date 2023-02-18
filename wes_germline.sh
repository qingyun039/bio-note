#!/bin/bash

function entryStep {
	step=$(( $step + 1 ))
	stepName=$1
	shift
	stepDir="$basedir/$step.$stepName"
	stepOk="$stepDir/ok"
	echo -e "[`date`] Processing $step.$stepName ..."
	echo -e "---------------[ `date`: $stepName ]---------------" >> $logfile
	if [[ ! -f $stepOk ]]
	then
		mkdir -p $stepDir
		oldDir=$PWD
		cd $stepDir
		{ $@ 2>> $logfile; }
		if [[ $? -eq 0 ]]
		then
			echo $tmpout > ok
		else
			exit 2
		fi
		cd $oldDir
	else
		tmpout=`cat $stepOk`
	fi
}

function fastqc {

	tmpout="$PWD/${sample}_1.fastq.gz $PWD/${sample}_2.fastq.gz"

	fastp -i $read1 -I $read2 -o ${sample}_1.fastq.gz -O ${sample}_2.fastq.gz \
        --unpaired1 ${sample}_unpair_1.fastq.gz \
        --unpaired2 ${sample}_unpair_2.fastq.gz \
        --failed_out ${sample}_failed.fastq.gz \
        -w $thread $fastp
}

function alignPlus {
	tmpin=$tmpout
        tmpout="$PWD/$sample.sorted.markdup.bam"

        bwa mem $bwa -R "@RG\tID:$sample\tSM:$sample\tPL:$platform" -t $thread $reference $tmpin |\
                samtools view -@ $thread -bS > $sample.bam && \
        sambamba sort -t $thread -o $sample.sorted.bam $sample.bam && \
        sambamba markdup -t $thread --overflow-list-size 600000 $sample.sorted.bam $sample.sorted.markdup.bam
}

function dataPreprocess {

	tmpin=$tmpout
	tmpout="$PWD/$sample.sorted.markdup.bqsr.bam"

	bwa mem -Y -R "@RG\tID:$sample\tSM:$sample\tPL:illumina" -t $thread $reference $tmpin |\
		samtools view -@ $thread -bS > $sample.bam && \
	sambamba sort -t $thread -o $sample.sorted.bam $sample.bam && \
	sambamba markdup -t $thread $sample.sorted.bam $sample.sorted.markdup.bam && \
	gatk BaseRecalibrator $interval_param \
		-R $reference \
		-I $sample.sorted.markdup.bam \
		-O $sample.sorted.markdup.bqsr.table \
		--known-sites $dbsnp \
		--known-sites $mills_and_1000g \
		--known-sites $phase1_indels && \
	gatk ApplyBQSR $interval_param \
		-R $reference \
		-I $sample.sorted.markdup.bam \
		-O $sample.sorted.markdup.bqsr.bam \
		--bqsr-recal-file $sample.sorted.markdup.bqsr.table
}

function parallelBQSR {
	tmpin=$tmpout
	tmpout="$PWD/$sample.bqsr.bam"

	dict=${reference%.*}.dict
	python /data/home/chenshulin/project/GERMLINE/software/gatk4-data-processing/split_group.py $dict && \
	perl -alne "print 'gatk BaseRecalibrator -R $reference -I $tmpin -O $sample.' . \$. . '.table --known-sites $dbsnp --known-sites $mills_and_1000g --known-sites $phase1_indels -L '. join(' -L ', @F)" sequence_grouping.txt | \
       	parallel --pipe -N1 bash && \
	gatk GatherBQSRReports \
		-I `ls $sample.*.table | sed ':a;N;$!ba;s/\n/ -I /g'` \
		-O $sample.table && \
	perl -alne "print 'gatk ApplyBQSR -R $reference -I $tmpin -O $sample.bqsr.'. \$. . '.bam -bqsr $sample.table -L '. join(' -L ', @F)" sequence_grouping.txt | \
	parallel --pipe -N1 bash && \
	gatk GatherBamFiles \
		-I `ls $sample.bqsr.*.bam | sed ':a;N;$!ba;s/\n/ -I /g'` \
		-O $sample.all.bam && \
	sambamba sort -t $thread -o $tmpout $sample.all.bam && sambamba index -t $thread $tmpout
}

function parallelHC {

	tmpin=$tmpout

	subfix='.vcf.gz'
	gvcf_param=''
	if [[ $gvcf -eq 1 ]]; then subfix=".g.vcf.gz"; gvcf_param="-ERC GVCF"; fi
	tmpout="$PWD/${sample}${subfix}"
	gatk SplitIntervals $interval_param \
		-R $reference \
		-scatter $thread \
		-O . && \
	ls *.interval_list | parallel gatk HaplotypeCaller $gvcf_param $hc \
		-R $reference \
		-I $tmpin \
		-L {} \
		-bamout $sample.{#}.bam \
		-O $sample.{#}${subfix} && \
	gatk MergeVcfs \
		-I `ls $sample.*${subfix} | sed ':a;N;$!ba;s/\n/ -I /g'` \
		-O $tmpout && \
	samtools merge -@ $thread -f $sample.bam $sample.*.bam && samtools index $sample.bam
}

# call
function callVariants {

	tmpin=$tmpout

	subfix='.vcf.gz'
	gvcf_param=''
	if [[ $gvcf -eq 1 ]]; then subfix=".g.vcf.gz"; gvcf_param="-ERC GVCF"; fi
	tmpout="$PWD/${sample}${subfix}"
	gatk HaplotypeCaller $gvcf_param $interval_param $hc \
		-R $reference \
		-I $tmpin \
		-O ${sample}${subfix} \
		-bamout $sample.hc.bam
}

function strelkaGermline {
	tmpin=$tmpout
	tmpout="$PWD/results/variants/variants.vcf.gz"
	bedtools sort -i $bedfile > region.bed && \
	bgzip -c region.bed > region.bed.gz && \
	tabix -p bed region.bed.gz && \
	/data/software/src/strelka-2.9.10.centos6_x86_64/bin/configureStrelkaGermlineWorkflow.py \
		--bam $tmpin \
		--referenceFasta $reference \
		--runDir . \
		--callRegions region.bed.gz \
		--exome && \
	./runWorkflow.py  -m local -j $thread
}

function varscanGermline {
	tmpin=$tmpout
	samtools mpileup -B -f $reference $tmpin | java -jar /data/home/chenshulin/.local/src/VarScan.v2.4.2.jar mpileup2snp --output-vcf 1 > varscan.snp.vcf && \
	samtools mpileup -B -f $reference $tmpin | java -jar /data/home/chenshulin/.local/src/VarScan.v2.4.2.jar mpileup2indel --output-vcf 1 > varscan.indel.vcf

}

#function jointtype {

#}

function filterCNN {

	tmpin=$tmpout
	tmpout="$PWD/${sample}.vcf"
	gatk CNNScoreVariants $interval_param \
		-R $reference \
		-V $tmpin \
		-O ${sample}_cnn_annotated.vcf.gz && \
	gatk FilterVariantTranches \
		-V ${sample}_cnn_annotated.vcf.gz \
		--output ${sample}.vcf \
		-resource ${hapmap} \
		-resource ${mills_and_1000g} \
		-resource ${phase1_snps} \
		-info-key CNN_1D \
		--snp-tranche 99.95 \
		--indel-tranche 99.4
}

function filterHard {
	tmpin=$tmpout
	tmpout="$PWD/$sample.vcf"
	gatk SelectVariants \
		-R $reference \
		-V $tmpin \
		-select-type SNP \
		-O $sample.snps.vcf.gz && \
	gatk VariantFiltration \
		-R $reference \
		-V $sample.snps.vcf.gz \
		-O $sample.snps.filtered.vcf.gz \
		--filter-expression "QD < 6.5" \
		--filter-name "QD_lt_6.5" \
		--filter-expression "QUAL < 30.0" \
		--filter-name "QUAL_lt_30" \
		--filter-expression "FS > 100.0" \
		--filter-name "FS_gt_100" \
		--filter-expression "MQ < 40.0" \
		--filter-name "MQ_lt_40" \
		--filter-expression "MQRankSum < -12.5" \
		--filter-name "MQRS_lt_n12.5" \
		--filter-expression "ReadPosRankSum < -8.0" \
		--filter-name "RPRS_lt_n8" && \
	gatk SelectVariants \
		-R $reference \
		-V $tmpin \
		-select-type INDEL \
		-O $sample.indels.vcf.gz && \
	gatk VariantFiltration \
		-R $reference \
		-V $sample.indels.vcf.gz \
		-O $sample.indels.filtered.vcf.gz \
		--filter-expression "QD < 4.0" \
		--filter-name "QD_lt_4" \
		--filter-expression "QUAL < 200.0" \
		--filter-name "QUAL_lt_200" \
		--filter-expression "FS > 50.0" \
		--filter-name "FS_gt_50" \
		--filter-expression "ReadPosRankSum < -5.0" \
		--filter-name "RPRS_lt_n5" && \
	gatk MergeVcfs \
		-I $sample.snps.filtered.vcf.gz \
		-I $sample.indels.filtered.vcf.gz \
		-O $sample.vcf
}

function annoAnnovar {
	tmpin=$tmpout
	tmpout="$PWD/$sample.${buildver}_multianno.txt"
	#vt decompose $tmpin | vt normalize -r $reference - > $sample.norm.vcf && \
	perl /data/database/annovar/table_annovar.pl $tmpin $humandb \
		-out $sample \
		-buildver ${buildver} \
		-remove \
		-protocol refGene,cytoBand,avsnp150,1000g2015aug_all,1000g2015aug_eas,clinvar,intervar,dbnsfp41a,gnomad211_exome \
		-operation g,r,f,f,f,f,f,f,f \
		-nastring . \
		-vcfinput -polish
}

function Test {

	echo job.{#}
}


if [[ $# -eq 1 ]];
then
	confile=$1
	source $confile
else
	echo "Read config from ENV..."
fi


#检查配置文件
if [[ -n $reference && -e $reference && -n $read1 && -e $read1 && -n $read2 && -e $read2 ]]
then
	echo -e "INPUT:\n\treference: $reference\n\tbedfile: $bedfile\n\tread1: $read1\n\tread2: $read2"
else
	echo "config must include: reference, read1, read2"
	exit 2
fi

guess=${read1##*/}
guess=${guess%%.*}
#其他默认的配置
echo "
Basic:
	sample: ${sample:=$guess}
	basedir: ${basedir:=.}
	platform: ${platform:=illumina}
	thread: ${thread:=36}
	gvcf: ${gvcf:=0}
	logfile: ${logfile:=$PWD/$sample.log}
Params:
	fastp: ${fastp:="-q 15"}
	bwa: ${bwa:="-Y"}
	hc: ${hc:=""}
Database:
	dbsnp: ${dbsnp:=/data/database/gatk-bundle-hg19/dbsnp_138.hg19.vcf}
	mills_and_1000g: ${mills_and_1000g:=/data/database/gatk-bundle-hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf}
	phase1_indels: ${phase1_indels:=/data/database/gatk-bundle-hg19/1000G_phase1.indels.hg19.sites.vcf}
	phase1_snps: ${phase1_snps:=/data/database/gatk-bundle-hg19/1000G_phase1.snps.high_confidence.hg19.sites.vcf}
	hapmap: ${hapmap:=/data/database/gatk-bundle-hg19/hapmap_3.3.hg19.sites.vcf}
	humandb: ${humandb:=/data/home/chenshulin/humandb/}
	buildver: ${buildver:=hg19}
"

step=0
interval_param=''
if [[ -n $bedfile ]]; then interval_param="-L $bedfile"; fi

source /data/home/chenshulin/.local/miniconda3/etc/profile.d/conda.sh

conda activate gatk

echo -e "---------------[ `date`: BEGIN ]---------------" > $logfile

entryStep fastqc fastqc;

entryStep align alignPlus;

entryStep bqsr parallelBQSR;

#entryStep varscan varscanGermline;

entryStep calls parallelHC;

if [[ $gvcf -ne 1 ]]; then

entryStep filter filterCNN;

entryStep hardf filterHard;

entryStep anno annoAnnovar;

fi

echo -e "---------------[ `date`: END ]---------------" >> $logfile
