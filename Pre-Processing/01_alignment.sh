for fastq in `cat $FASTQDIRS`; do
  id=`echo $fastq | cut -d/ -f9 | cut -d\. -f1 | cut -d_ -f3-15`
  echo $id >> $IDLIST
  ln -s $fastq ${FQSYM}
done

cat $IDLIST | sort | uniq > ${IDLIST}.tmp;
mv ${IDLIST}.tmp ${IDLIST};


for id in `cat $IDLIST`; do
  echo $id
  FQIDSYM=${FQSYM}/${id}/
  mkdir -p ${FQIDSYM}
  LIBRARY="LB:$LIBID"
  LANE="PU:L001"
  PLATFORM="PL:Illumina"
  SAMPLE="SM:$id"
  sID="ID:$id"
  ID=$id
  FASTQ1="${id}_L001_R1_001.fastq.gz"
  FASTQ2="${id}_L001_R2_001.fastq.gz"
  RGID="\"@RG\\t${sID}\\t${SAMPLE}\\t${PLATFORM}\\t${LANE}\\t${LIBRARY}\""
  # RGID=`cat RGID.txt`
  if [[ ! -f $FASTQ1 && ! -f $FASTQ2 ]]; then
    `zless ${FQSYM}/*${id}*R1*.fastq.gz | gzip > ${FQIDSYM}/${FASTQ1}`
    `zless ${FQSYM}/*${id}*R2*.fastq.gz | gzip > ${FQIDSYM}/${FASTQ2}`
  fi
  echo "bwa mem -M -t4 -R ${RGID} ${BWAINDEX} ${FQIDSYM}/${FASTQ1} ${FQIDSYM}/${FASTQ2} > ${BAMDIR}/${ID}.sam \\
  2> ${BAMDIR}/${ID}_bwamem.err" > $id.bwa.sh
  chmod +x $id.bwa.sh
  if [[ ! -f ${BAMDIR}/${ID}.bam ]]; then
    echo "$SAMTOOLS view -bhS ${BAMDIR}/${ID}.sam | samtools sort -@4 - ${BAMDIR}/${ID}.bam; rm ${BAMDIR}/${ID}.sam" >> $id.bwa.sh
  fi
  `source $id.bwa.sh`
done
