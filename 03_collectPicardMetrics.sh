


for id in `cat ${IDLIST}`; do
  PCFL=$PCQCDIR/${id}
  # mkdir -p $PCDIR
  BAM="${BAMDIR}/${id}.bam.bam"
  echo "$java -Xmx6g -jar $picard_jar CollectWgsMetrics \\
    I=$BAM.cocleaned.bam \\
    O=$PCFL.wgsMetrics.txt \\
    R=$REF" > ${id}.picard_metrics.sh
    echo "$java -Xmx6g -jar $picard_jar CollectInsertSizeMetrics \\
    I=$BAM.cocleaned.bam \\
    O=$PCFL.isize.txt \\
    H=$PCFL.histogram.pdf \\
    M=0.5" >> ${id}.picard_metrics.sh
    chmod +x ${id}.picard_metrics.sh
done
source ${id}.picard_metrics.sh
