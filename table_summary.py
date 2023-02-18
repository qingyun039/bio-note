#!/usr/bin/env python
import sys
import os
import json
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from collections import OrderedDict

def long_table(filename, skip = 0, delimiter = None):
    with open(filename) as f:
        linenum = 0
        rd = OrderedDict()
        for line in f:
            linenum += 1
            if linenum <= skip:
                continue
            line = line.split(delimiter)
            rd[line[0].strip()] = line[1].strip()
    return rd

def wide_table(filename, skip = 0, delimiter = None):
    with open(filename) as f:
        linenum = 0
        rd = OrderedDict()
        keys = list()
        values = list()
        for line in f:
            linenum += 1
            if linenum <= skip:
                continue
            line = line.split(delimiter)
            if keys:
                values = line
                break
            else:
                keys = line
    rd = OrderedDict(zip(keys, values))
    return rd

def fastp_json(filename):
    with open(filename) as f:
        data = json.loads(f.read())
    rd = OrderedDict()
    for k, v in data['summary']['before_filtering'].items():
        rd[k] = v
    for k in ["low_quality_reads", "too_many_N_reads", "too_short_reads"]:
        rd[k] = data["filtering_result"][k]
        rd[k+'_rate'] = data["filtering_result"][k] / data['summary']['before_filtering']['total_reads']
    for k in ["adapter_trimmed_reads", "adapter_trimmed_bases"]:
        rd[k] = data["adapter_cutting"][k]
        if k == "adapter_trimmed_reads":
            rd[k+'_rate'] = data["adapter_cutting"][k] / data['summary']['before_filtering']['total_reads']
        else:
            rd[k+'_rate'] = data["adapter_cutting"][k] / data['summary']['before_filtering']['total_bases']
    for k, v in data['summary']['after_filtering'].items():
        rd[k+'_clean'] = v
    return rd

def mgiseq_qc(fqdir, lane, barcode):
    fqstats = glob.glob(fqdir +'/' + lane + '/*_' + barcode + '_*.fq.fqStat.txt')
    print(barcode, lane)
    print(fqstats)
    l = len(fqstats)
    rd = OrderedDict()
    record = ['ReadNum', 'BaseNum', 'N_Count', 'N_Count%', 'GC%', 'Q20%', 'Q30%', 'EstErr%']
    for k in record:
        rd[k] = 0
    for filename in fqstats:
        with open(filename) as f:
            for line in f:
                line = line.strip().split()
                for k in record:
                    if line[0] == '#'+k:
                        rd[k] += float(line[1])
                        if k == 'N_Count':
                            rd['N_Count%'] = float(line[2])
    for k in record:
        if k[-1] == '%':
            rd[k] = rd[k] / l
    return(rd)


#
def bamqc(bamqc_dir, fqdir='/data/RawData/MGISEQ-2000/F300005191/L02/', lane = 'L01', barcode='DY-289'):
    bamdst_coverage = os.path.join(bamqc_dir, 'coverage.report')
    cov_summary_file = os.path.join(bamqc_dir, 'cov_summary.txt')
    rmdup_cov_summary_file = os.path.join(bamqc_dir, 'cov_summary_rmdup.txt')
    gc_summary_file = os.path.join(bamqc_dir, 'gc_bias.summary_metrics')
    hs_metrics_file = os.path.join(bamqc_dir, 'hybird_capture.metrics')

    result = mgiseq_qc(fqdir, lane, barcode)

    tmprs = fastp_json(os.path.join(bamqc_dir,'../../1.fastqc/fastp.json'))
    result.update(tmprs)

    tmprs = long_table(bamdst_coverage, 3, '\t')
    result.update(tmprs)

    tmprs = long_table(cov_summary_file, delimiter = '\t')
    result.update(tmprs)

    tmprs = long_table(rmdup_cov_summary_file, delimiter = '\t')
    tmprs =  OrderedDict([ (k+'rmdup', v) for k,v in tmprs.items() ])
    result.update(tmprs)

    tmprs = wide_table(gc_summary_file, 6, delimiter = '\t')
    result.update(tmprs)

    tmprs = wide_table(hs_metrics_file, 6, delimiter = '\t')
    result.update(tmprs)

    return(result)

def qc_plot(bamqc_dirs, filename, head=375):
    dat_files = [os.path.join(i, filename) for i in bamqc_dirs]
    #df = pd.DataFrame()
    fig = plt.figure()
    ax = fig.add_axes([0,0,1,1])
    ax2 = ax.twinx()
    for f in dat_files:
        dat = pd.read_table(f, header=None)
        #df = pd.concat([df, dat.iloc[:head, :]])
        dat = dat.iloc[:head, :]
        sns.lineplot(data=dat, x=0, y=2, ax=ax, color='g')
        sns.lineplot(data=dat, x=0, y=4, ax=ax2, color='b')
    ax.set_title(filename, fontweight='bold')
    ax.set_xticks(list(range(0, head, 50)))
    ax.set_xlabel('depth')
    ax.set_ylim(0, 0.01)
    ax.set_ylabel('freq')
    ax.figure.savefig('test2.png')

############## main #############
if __name__ == '__main__':
    try:
        ana_dir = sys.argv[1]
        fqdir = sys.argv[2]
    except:
        print("Usage: {} <Analysis DIR>".format(sys.argv[0]))

    bamqc_dirs = glob.glob(os.path.join(ana_dir, '*', '2.align', 'bamqc'))
    bamqc_dirs = sorted(bamqc_dirs)
    samples = [x.split('/')[-3] for x in bamqc_dirs]
    samples = [ sample.split('_') for sample in samples]
    samples = pd.DataFrame(samples)
    samples.columns = ['sample', 'barcode', 'hs']

    all_aggr = list()
    for bamqc_dir in bamqc_dirs:
        barcode, lane = bamqc_dir.split('/')[-3].split('_')[1:3]
        lane = lane[0:3]
        if True:
            all_aggr.append(bamqc(bamqc_dir, fqdir, lane, barcode))
        else:
            all_aggr.append(bamqc(bamqc_dir))

    df = pd.DataFrame(all_aggr)
    writer = pd.ExcelWriter('qc_table.xlsx')
    cols = [0, 8, 13, 24, 57, 63, 98, 101, 104]
    pd.concat([samples, df.iloc[:,cols]], axis=1).to_excel(writer, 'Some Data')
    pd.concat([samples, df], axis=1).to_excel(writer, 'All Data')


    #cols = [0,1,2,3,4,5,6,7,8,9,12,13,23,24,25,26,27,28,31,32,36,41,57,61,63,67,70,71,91,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,173]
    writer.save()

    # bar plot
    cols = ['[Target] Fraction of Target Reads in all reads', '[Target] Fraction of Target Data in all data', 'Uniformity(0.2xAvg %)', 'Uniformity(0.5xAvg %)',
        'Coverage >=1x(%)', 'Coverage >=5x(%)', 'Coverage >=10x(%)', 'Coverage >=20x(%)', 'Coverage >=30x(%)']
    ncols = ['On-Target', 'On-Target Base', '.2Uniform', '.5Uniform', '>=1X', '>=5X', '>=10X', '>=20X', '>=30X']
    #cols = ['total_reads','total_reads_clean', '[Total] Mapped Reads', '[Total] MapQuality above cutoff reads']
    #ncols = ['total_reads', 'clean_reads', 'mapped reads', 'mapq>20 reads']
    #cols = ['low_quality_reads_rate', 'too_many_N_reads_rate', 'adapter_trimmed_reads_rate']
    #ncols = ['qulity', 'N', 'adapter']
    #cols = ['adapter_trimmed_reads_rate']
    #ncols = ['adapter']
    #cols = ['q20_rate', 'q30_rate', 'q20_rate_clean', 'q30_rate_clean']
    #ncols = ['q20', 'q30', 'q20_clean', 'q30_clean']
    #cols = ['[Total] Fraction of PCR duplicate reads']
    #ncols = ['duplicate rate']
    bardf = df[cols]
    bardf['[Target] Fraction of Target Reads in all reads'] = bardf['[Target] Fraction of Target Reads in all reads'].str.replace('%', '')
    bardf['[Target] Fraction of Target Data in all data'] = bardf['[Target] Fraction of Target Data in all data'].str.replace('%', '')
    #bardf['[Total] Fraction of PCR duplicate reads'] = bardf['[Total] Fraction of PCR duplicate reads'].str.replace('%', '')
    bardf = bardf.apply(pd.to_numeric)
    bardf.rename(columns = dict(zip(cols, ncols)), inplace = True)
    bardf['sample'] = samples['sample']
    bardf['barcode'] = samples['barcode'].apply(lambda x: x.split('-')[0])
    bardf['hs'] = samples['hs'].apply(lambda x: 'NAD' if x == 'L04' else 'IGT')
    bardf['group'] = bardf[['barcode', 'hs']].apply(lambda x: '-'.join(x), axis=1)
    bardf = pd.melt(bardf, id_vars = ['sample', 'group'], value_vars=ncols)
    #print(bardf)
    sns.set_theme(style="whitegrid")
    ax = sns.barplot(x='variable', y = 'value', hue = 'group', data = bardf)
    ax.set_xticklabels(ncols, rotation=15)
    ax.figure.savefig('ontarget.png')

    #qc_plot(bamqc_dirs, 'depth_distribution.plot')


