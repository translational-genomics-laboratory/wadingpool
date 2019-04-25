# module load picard/1.90
module load samtools/1.2
# module load igenome-human/hg19
module load bwa/0.7.9a
module load gatk/3.5.0

for id in `cat $IDLIST`; do
  BAM="${BAMDIR}/${id}.bam.bam"
  if [[ -f $BAM ]]; then
    echo "$java -Xmx4g -jar $picard_jar MarkDuplicates AS=TRUE \\
      I=$BAM \\
      M=$BAM.dup_metric.txt \\
      O=$BAM.dup.bam \\
      2> $LOGDIR/markduplicates.log" > ${id}.coclean.sh

    echo "$SAMTOOLS index $BAM.dup.bam 2>> $LOGDIR/indelrealign.log" >> ${id}.coclean.sh
    echo "$java -Xmx4g -Djava.io.tmpdir=tmpdir/ \\
    -jar $gatk_dir/GenomeAnalysisTK.jar \\
    -T RealignerTargetCreator -nt 8 \\
      -I $BAM.dup.bam \\
      -R $REF \\
      -o $BAM.intervals \\
      -dt NONE \
      -known $grch371000gIndels \\
      -known $grch37MillsIndels \\
      2>> $LOGDIR/indelrealign.log" >> ${id}.coclean.sh

      echo "$java -Xmx4g -Djava.io.tmpdir=tmpdir/ \\
      -jar $gatk_dir/GenomeAnalysisTK.jar \\
      -T IndelRealigner \
      -I $BAM.dup.bam \\
      -R $REF \\
      -targetIntervals $BAM.intervals \\
      -dt NONE \\
      -o $BAM.dup.realign.bam \\
      -known $grch371000gIndels \\
      -known $grch37MillsIndels \\
      2>> $LOGDIR/indelrealign.log" >>${id}.coclean.sh

     echo "$java -Xmx4g -Djava.io.tmpdir=tmpdir/ \\
     -jar $gatk_dir/GenomeAnalysisTK.jar \\
     -T BaseRecalibrator -nct 8 \\
      -I $BAM.dup.realign.bam \\
      -o $BAM.dup.realign.recal.grp \\
      -R $REF \\
      -knownSites $dbsnpVcf \\
      -rf BadCigar \\
      -cov ReadGroupCovariate \\
      -cov ContextCovariate \\
      -cov CycleCovariate \\
      -cov QualityScoreCovariate \\
      -dt NONE \\
      -knownSites $grch371000gIndels \\
      -knownSites $grch37MillsIndels \\
      2>> $LOGDIR/recalibrate.log" >> ${id}.coclean.sh

    echo "$java -Xmx4g -Djava.io.tmpdir=tmpdir/ \\
    -jar $gatk_dir/GenomeAnalysisTK.jar \\
    -T PrintReads -nct 8 \\
      -R $REF \\
      -I $BAM.dup.realign.bam \\
      -BQSR $BAM.dup.realign.recal.grp \
      -o $BAM.cocleaned.bam \\
      -rf BadCigar -dt NONE \\
      2>> $LOGDIR/recalibrate.log">> ${id}.coclean.sh

      echo "rm $BAM.dup.bam ;
      rm $BAM.dup.bam.bai ;
      rm $BAM.intervals ;
      rm $BAM.dup.realign.bam ;
      rm $BAM.dup.realign.bai ;
      rm $BAM.dup.realign.recal.grp ;" >> ${id}.coclean.sh
  fi
  chmod +x ${id}.coclean.sh
done
source ${id}.coclean.sh
