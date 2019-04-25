for id in `cat ${IDLIST}`; do
  PCFL=$PCQCDIR/${id}
  # mkdir -p $PCDIR
  BAM="${BAMDIR}/${id}.bam.bam"
  echo "telomerecat bam2telbam \\
  -v 2 \\
  -p 4 \\
  $BAM.cocleaned.bam" > ${id}.telomereBam.sh
  echo "mv ${id}.bam.bam.cocleaned_telbam.bam $TLDIR/${id}.telbam.bam">> ${id}.telomereBam.sh
  chmod +x ${id}.telomereBam.sh
  source ${id}.telomereBam.sh
done


echo "telomerecat telbam2length \\
-p 4 \\
-v 2 \\
-N 100 \\
--output $TL/allTelbamLengths.csv \\
\$(ls $TLDIR/*telbam.bam)" > telomereLength.sh

chmod +x telomereLength.sh
source telomereLength.sh
