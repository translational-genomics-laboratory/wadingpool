## symlink FASTQ directory and align
bold=$(tput bold)
normal=$(tput sgr0)
cat << EOF
${bold}WadingPool:${normal} Coclean bam files
    Helper script to take in the bwa aligned and demultiplexed bam files and run them through GATK post-processing pipeline
    Assumes the fastq.gz files take the following identification structure:
      
      CHX-022-T_S44_L001_R1_001.fastq.gz
      [ID]_[LANE]_[READ]_[MULTIPLEX_ID].fastq.gz
      ID:CHX-022-T_S44
EOF


mkdir -p ${PDIR}/data/sh_scripts/coclean
cd ${PDIR}/data/sh_scripts/coclean

rm id_list.txt cocleanBam.sh
ln -s ${IDLIST} .
ln -s ${GIT}/preprocessing/alignment/coclean/cocleanBam.sh .

echo "${bold}WadingPool:${normal} cocleaning bam files...";
cat << EOF > queueJobs.sh
for i in \$(cat id_list.txt); do
  qsub cocleanBam.sh \$i ${PDIR} ${REFERENCEFILES}
done
EOF
sh queueJobs.sh
