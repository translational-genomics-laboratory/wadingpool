#!/bin/bash
#
#$ -cwd
#$ -S /bin/bash
#
############################################################
bold=$(tput bold)
normal=$(tput sgr0)

## Specify PDIR and LIBID as passed in or hardcoded
if [ $# -eq 0 ]; then
  echo "${bold}WadingPool:${normal} Using hardcoded PDIR and REFERENCEFILES"
  BAM='CHX_test'
  PDIR='/mnt/work1/users/pughlab/projects/test_dir/swgs'
  REFERENCEFILES='/mnt/work1/users/pughlab/references/igenome/coclean_reference_files.sh'
else
  echo "${bold}WadingPool:${normal} Using PDIR and REFERENCEFILES specified by configuration file"
  BAM=$1
  PDIR=$2
  REFERENCEFILES=$3
fi
source ${REFERENCEFILES}
BAM=${BAM}.bam


FASTQDIR=$PDIR'/data/fastq'
NBAMDIR=$PDIR'/data/bam'
CBAMDIR=$PDIR'/data/cocleaned_bam'

module load picard/1.9.1
module load samtools/1.2
module load igenome-human/hg19
module load bwa/0.7.9a
module load gatk/3.5

##############################
## M1: Mark duplicates on the realigned BAM file
stripBAM=${BAM//.bam/}
LOGDIR=$PDIR'/data/logs/coclean/'$stripBAM
mkdir -p $LOGDIR
CBAMDIR=$CBAMDIR'/'$stripBAM
mkdir -p $CBAMDIR

echo "["$(date)"] M1: Marking duplicates for "$BAM > $LOGDIR'/reprocess.out'
if [ -f $NBAMDIR'/'$BAM ]; then
    java -Xmx4g -jar $picard_dir/MarkDuplicates.jar AS=TRUE \
      I=$NBAMDIR'/'$BAM \
      METRICS_FILE=$CBAMDIR'/'$stripBAM'.dup_metric.txt' \
      O=$CBAMDIR'/'$stripBAM'.dup.bam' \
      2> $LOGDIR'/markduplicates.log'
else
    echo "Error:  file not present - "$NBAMDIR'/'$BAM  > $LOGDIR'/reprocess.err'
fi

##############################
## M2: Indel realignment
echo "["$(date)"] M2: Realigning for indels" >> $LOGDIR'/reprocess.out'
if [ -f $CBAMDIR'/'$stripBAM'.dup.bam' ]; then
    echo -e "\t Indexing duplicated marked BAM file.." >> $LOGDIR'/reprocess.out'
    samtools index $CBAMDIR'/'$stripBAM'.dup.bam' 2>> $LOGDIR'/indelrealign.log'

    echo -e "\t Creating the interval targets for indel realignment.." >> $LOGDIR'/reprocess.out'
    java -Xmx4g -Djava.io.tmpdir=tmpdir/ -jar $gatk_dir/GenomeAnalysisTK.jar -T RealignerTargetCreator -nt 8 \
      -I $CBAMDIR'/'$stripBAM'.dup.bam' \
      -R $REF \
      -o $CBAMDIR'/'$stripBAM'.intervals' \
      -dt NONE \
      -known $grch371000gIndels \
      -known $grch37MillsIndels \
      2>> $LOGDIR'/indelrealign.log'
else
    echo "Error:  file not present - "$CBAMDIR'/'$stripBAM'.dup.bam'  >> $LOGDIR'/reprocess.err'
fi

if [ -f $CBAMDIR'/'$stripBAM'.intervals' ]; then
    echo -e "\t Realigning according to the "$CBAMDIR'/'$stripBAM'.intervals'" targets..." >> $LOGDIR'/reprocess.out'
    java -Xmx4g -Djava.io.tmpdir=tmpdir/ -jar $gatk_dir/GenomeAnalysisTK.jar -T IndelRealigner \
      -I $CBAMDIR'/'$stripBAM'.dup.bam' \
      -R $REF \
      -targetIntervals $CBAMDIR'/'$stripBAM'.intervals' \
      -dt NONE \
      -o $CBAMDIR'/'$stripBAM'.dup.realign.bam' \
      -known $grch371000gIndels \
      -known $grch37MillsIndels \
      2>> $LOGDIR'/indelrealign.log'
else
    echo "Error:  file not present - "$CBAMDIR'/'$stripBAM'.intervals'  >> $LOGDIR'/reprocess.err'
fi

##############################
## M3: Base Recalibration
echo "["$(date)"] M3: Recalibrating bases ..." >> $LOGDIR'/reprocess.out'
if [ -f $CBAMDIR'/'$stripBAM'.dup.realign.bam' ]; then    
    echo -e "\t Calculating the recalibrated bases..." >> $LOGDIR'/reprocess.out'
    java -Xmx4g -Djava.io.tmpdir=tmpdir/ -jar $gatk_dir/GenomeAnalysisTK.jar -T BaseRecalibrator -nct 8 \
      -I $CBAMDIR'/'$stripBAM'.dup.realign.bam' \
      -o $CBAMDIR'/'$stripBAM'.dup.realign.recal.grp' \
      -R $REF \
      -knownSites $dbsnpVcf \
      -rf BadCigar \
      -cov ReadGroupCovariate \
      -cov ContextCovariate \
      -cov CycleCovariate \
      -cov QualityScoreCovariate \
      -dt NONE \
      -knownSites $grch371000gIndels \
      -knownSites $grch37MillsIndels \
      2>> $LOGDIR'/reacalibrate.log'
else
    echo "Error:  file not present - "$CBAMDIR'/'$stripBAM'.dup.realign.bam'  >> $LOGDIR'/reprocess.err'
fi

if [ -f $CBAMDIR'/'$stripBAM'.dup.realign.recal.grp' ]; then    
    echo -e "\t Printing recalibrated bam..." >> $LOGDIR'/reprocess.out'
    java -Xmx4g -Djava.io.tmpdir=tmpdir/ -jar $gatk_dir/GenomeAnalysisTK.jar -T PrintReads -nct 8 \
      -R $REF \
      -I $CBAMDIR'/'$stripBAM'.dup.realign.bam' \
      -BQSR $CBAMDIR'/'$stripBAM'.dup.realign.recal.grp' \
      -o $CBAMDIR'/'$stripBAM'.cocleaned.bam' \
      -rf BadCigar -dt NONE \
      2>> $LOGDIR'/reacalibrate.log'
else
    echo "Error:  file not present - "$CBAMDIR'/'$stripBAM'.dup.realign.recal.grp'  >> $LOGDIR'/reprocess.err'
fi

##############################
## M4: Cleanup
echo "["$(date)"] M4: Cleaning up intermediate files..." >> $LOGDIR'/reprocess.out'
if [ -f $CBAMDIR'/'$stripBAM'.cocleaned.bam' ]; then    
    rm $CBAMDIR'/'$stripBAM'.dup.bam'
    rm $CBAMDIR'/'$stripBAM'.dup.bam.bai'
    rm $CBAMDIR'/'$stripBAM'.intervals'
    rm $CBAMDIR'/'$stripBAM'.dup.realign.bam'
    rm $CBAMDIR'/'$stripBAM'.dup.realign.bai'
    rm $CBAMDIR'/'$stripBAM'.dup.realign.recal.grp'
    echo "["$(date)"] Completed cocleaning bams." >> $LOGDIR'/reprocess.out'
else
    echo "["$(date)"] Failed.  Did not remove intermediate files." >> $LOGDIR'/reprocess.out'
fi

