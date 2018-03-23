bold=$(tput bold)
normal=$(tput sgr0)

## Specify PDIR and LIBID as passed in or hardcoded
if [ $# -eq 0 ]; then
  echo "${bold}WadingPool:${normal} Using hardcoded PDIR and LIBID"
  PDIR='/mnt/work1/users/pughlab/projects/test_dir/swgs'
  LIBID='171120_NB501085_0175_AHNJTFAFXX_Milosevic_Rene-Pugh'
else
  echo "${bold}WadingPool:${normal} Using PDIR and LIBID specified by configuration file"
  PDIR=$1
  LIBID=$2
fi
PDIR=${PDIR}'/data'


## Generate shell scripts to merge and align all fastq files
rm ${PDIR}/sh_scripts/bwa/queueJobs.sh
for id in $(cat id_list.txt); do
  cat << EOF > bwa.$id.sh
#!/bin/bash
#
#$ -cwd
#$ -S /bin/bash
#

module load bwa/0.7.9a
module load samtools/1.2
module load igenome-human/hg19

# Read group info:
LIBRARY='LB:$LIBID'
LANE='PU:L001'
PLATFORM='PL:Illumina'
SAMPLE='SM:$id'
sID='ID:$id'

BWAINDEX='/mnt/work1/data/genomes/human/hg19/iGenomes/Sequence/BWAIndex/genome.fa'

ID=$id
FASTQ1="${id}_R1_001.fastq.gz"
FASTQ2="${id}_R2_001.fastq.gz"
RGID="@RG\\t\$sID\\t\$SAMPLE\\t\$PLATFORM\\t\$LANE\\t\$LIBRARY"


echo "Concatenating fastq files"
zcat \$(ls -1 ${PDIR}/fastq_symlinks/* | grep "${id}.*R1.*fastq.gz") | gzip > ${PDIR}/fastq/${id}_R1_001.fastq.gz
zcat \$(ls -1 ${PDIR}/fastq_symlinks/* | grep "${id}.*R2.*fastq.gz") | gzip > ${PDIR}/fastq/${id}_R2_001.fastq.gz

bwa mem -M -t4 \\
   -R "\${RGID}" \\
   \${BWAINDEX} \\
   ${PDIR}/fastq/\${FASTQ1} \\
   ${PDIR}/fastq/\${FASTQ2} \\
   > ${PDIR}/bam/\${ID}.sam \\
   2> ${PDIR}/bam/\${ID}_bwamem.err

if [ -f ${PDIR}/bam/\${ID}.sam ]; then
    samtools view -bhS ${PDIR}/bam/\${ID}.sam |\\
      samtools sort -@4 - ${PDIR}/bam/\${ID}
    
    samtools index ${PDIR}/bam/\${ID}.bam
else
    echo "Error:  file not present "
fi

if [ -f ${PDIR}/bam/\${ID}.sam ]; then
   rm ${PDIR}/bam/\${ID}.sam
   echo "["\$(date)"] Completed."
else
   echo "["\$(date)"] Failed.  Did not remove intermediate files."
fi
EOF

  echo "qsub bwa.${id}.sh" >> ${PDIR}/sh_scripts/bwa/queueJobs.sh
done

