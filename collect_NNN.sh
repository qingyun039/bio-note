#/bin/bash
# 统计测序中出现NNNNN的read的最高的芯片格子
for sequencer in $@;
do
    for run in `ls -d $sequencer/*`;
    do
        chip=${run##*/}
        echo "$run"
        for read in `ls $run/*/*.fq.gz | parallel -j 16 seqkit grep -p 'NNNNN' -s {} | grep "^@$chip"`
        do
            echo ${read:0:21}
        done | sort | uniq -c | sort -n -r | head -1
    done
done

