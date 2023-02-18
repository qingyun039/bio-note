
# boxplot
samtools view  320221459200/2.align/320221459200.sorted.markdup.bam | cut -f5 | gnuplot -p -e "set datafile separator '\t'; binwidth=5; bin(x,width)=width*floor(x/width);plot '-' using (bin(\$1,binwidth)):(1.0) smooth freq with boxes"

# depth distribution plot
head -n 375 320222146400_L3/2.align/qc/depth_distribution.plot | gnuplot -p -e "set datafile separator '\t'; plot '-' using 1:4 w l"

# insert size distribution plot
cat insertsize.plot | gnuplot -p -e "set datafile separator '\t'; plot '-' using 1:3 w l"

# region depth z-score
 zcat region.tsv.gz | perl -MStatistics::Basic=:all -MScalar::Util=looks_like_number -anle 'push @a, $F[3] if looks_like_number $F[3]; END{ $mean =  mean(@a); $sd = stddev(@a); print ++$index . "\t" . ($_ - $mean)/$sd for (@a)}' | gnuplot -p -e "set datafile separator '\t'; plot '-' using 1:2 w p"
