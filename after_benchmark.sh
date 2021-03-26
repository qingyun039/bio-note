#!/bin/bash

# 使用hap.py进行分析后，通过以下步骤分析假阳性，假阴性：
TRUTH=/data/home/chenshulin/project/GERMLINE/compare/truthSet/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf.gz
HIGHCONF=/data/home/chenshulin/project/GERMLINE/compare/truthSet/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel.bed
REF=/data/database/gatk-bundle-b37/human_g1k_v37_decoy.fasta

input=/data/home/chenshulin/project/GERMLINE/Analysis/MGIdemo/5.hardf/MGIdemo.vcf.gz
bedfile=/data/home/chenshulin/project/GERMLINE/Data/demo/MGI_Exome_Capture_V5_b37.bed
bamfile=/data/home/chenshulin/project/GERMLINE/Analysis/MGIdemo/3.bqsr/MGIdemo.bqsr.bam


# 使用python将bam-readcount的输出转换成XLSX格式
read -r -d '' pystr <<PYTHON
import sys
import xlsxwriter
posinfo = list()
for line in sys.stdin:
    line = line.strip().split()
    if len(line) <= 5:
        line = line + ['', '', '', '0',  '=:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00']
    k = line[0:4]
    k.append(line[7])
    info = list()
    for item in line[8:]:
        nuc = item.split(':')
        if int(nuc[1]) > 0:
            info.append(nuc)
    info = sorted(info, key=lambda x: int(x[1]), reverse=True)
    info.insert(0,k)
    posinfo.append(info)
posinfo = sorted(posinfo, key=lambda x: int(x[0][-1]))

workbook = xlsxwriter.Workbook(sys.argv[1])
worksheet = workbook.add_worksheet()
head_format = workbook.add_format({'bold': True, 'fg_color': '#D7E4BC',})

merge_format = workbook.add_format({
    'bold': 1,
    'border': 2,
    'align': 'center',
    'valign': 'vcenter',})
header = ['Chrom', 'Pos', 'Ref', 'Alt', 'Depth', 'base', 'count', 'mapq', 'qual', 'seq', 'plus', 'minus', 'pos_frac', 'mis_frac', 'mis_qual', 'q2', 'q2_reads', 'clip_len', '3p']

row = 0
col = 0
for i, s in enumerate(header):
    worksheet.write(row, col + i, s, head_format)

row = 1
for item in posinfo:
    n = len(item) - 1
    baseinfo = item[0]
#    print(row, row + n - 1)
    for i,s in enumerate(baseinfo):
        if n > 1:
            worksheet.merge_range(row, col + i, row + n - 1, col + i, s, merge_format)
        else:
            worksheet.write(row, col + i, s, merge_format)
    for i, info in enumerate(item[1:]):
        for j, s in enumerate(info):
            worksheet.write(row+i, col+j+5, s)
    if n == 0:
        row += 1
    else:
        row += n 
    col = 0
workbook.close()
PYTHON

bam_readcount(){
	REF=$1
	bamfile=$2
	shift 2
	echo $@ `bam-readcount -w 0 -f $REF $bamfile $1:$2-$2 2> /dev/null`
}
export -f bam_readcount
#bam_readcount 1 17009 A C

#hap.py $TRUTH  $input -o benchmark -r $REF -T $bedfile -f $HIGHCONF

#bcftools view -i 'TYPE="SNP" && BD[0] = "TP"' -o benchmark.TP.snp.vcf.gz -O z  benchmark.vcf.gz
#bcftools view -i 'TYPE="SNP" && BD[1] = "FP"' -o benchmark.FP.snp.vcf.gz -O z  benchmark.vcf.gz
#bcftools view -i 'TYPE="SNP" && BD[0] = "FN"' -o benchmark.FN.snp.vcf.gz -O z benchmark.vcf.gz
#bcftools view -i 'TYPE="INDEL" && BD[0] = "TP"' -o benchmark.TP.indel.vcf.gz -O z  benchmark.vcf.gz
#bcftools view -i 'TYPE="INDEL" && BD[1] = "FP"' -o benchmark.FP.indel.vcf.gz -O z benchmark.vcf.gz
#bcftools view -i 'TYPE="INDEL" && BD[0] = "FN"' -o benchmark.FN.indel.vcf.gz -O z benchmark.vcf.gz

#bcftools view -T benchmark.TP.snp.vcf.gz $input | bcftools query -f "%CHROM-%POS\t%REF\t%ALT\t%DP\t%QUAL\t%QD\t%FS\t%MQ\t%MQRankSum\t%ReadPosRankSum\t%SOR\tTP\tSNP\n" > benchmark.tsv
#bcftools view -T benchmark.FP.snp.vcf.gz $input | bcftools query -f "%CHROM-%POS\t%REF\t%ALT\t%DP\t%QUAL\t%QD\t%FS\t%MQ\t%MQRankSum\t%ReadPosRankSum\t%SOR\tFP\tSNP\n" >> benchmark.tsv
#bcftools view -T benchmark.FN.snp.vcf.gz $input | bcftools query -f "%CHROM-%POS\t%REF\t%ALT\t%DP\t%QUAL\t%QD\t%FS\t%MQ\t%MQRankSum\t%ReadPosRankSum\t%SOR\tFN_GT\tSNP\n" >> benchmark.tsv
bcftools query -i 'BD[0]="FN" && BD[1]="."' -f "%CHROM %POS %REF %ALT\n" benchmark.FN.snp.vcf.gz | \
	sed "s#^#$REF $bamfile #" | \
        parallel -a - --colsep " " bam_readcount | \
        python -c "$pystr" SNP_FN_NC.xlsx

#bcftools view -T benchmark.TP.indel.vcf.gz $input | bcftools query -f "%CHROM-%POS\t%REF\t%ALT\t%DP\t%QUAL\t%QD\t%FS\t%MQ\t%MQRankSum\t%ReadPosRankSum\t%SOR\tTP\tINDEL\n" >> benchmark.tsv
#bcftools view -T benchmark.FP.indel.vcf.gz $input | bcftools query -f "%CHROM-%POS\t%REF\t%ALT\t%DP\t%QUAL\t%QD\t%FS\t%MQ\t%MQRankSum\t%ReadPosRankSum\t%SOR\tFP\tINDEL\n" >> benchmark.tsv
#bcftools view -T benchmark.FP.indel.vcf.gz $input | bcftools query -f "%CHROM-%POS\t%REF\t%ALT\t%DP\t%QUAL\t%QD\t%FS\t%MQ\t%MQRankSum\t%ReadPosRankSum\t%SOR\tFN_GT\tINDEL\n" >> benchmark.tsv

bcftools query -i 'BD[0]="FN" && BD[1]="."' -f "%CHROM %POS %REF %ALT\n" benchmark.FN.indel.vcf.gz | \
	sed "s#^#$REF $bamfile #" | \
	parallel -a - --colsep " " bam_readcount | \
	python -c "$pystr" INDEL_FN_NC.xlsx
