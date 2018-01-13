## symlink FASTQ directory and align
bold=$(tput bold)
normal=$(tput sgr0)
cat << EOF
${bold}WadingPool:${normal} Telomerecat
    Helper script to take in all the non-cocleaned bam files and run telomerecat to quantify the length of the telomeres
    Assumes the fastq.gz files take the following identification structure:
      
      CHX-022-T_S44_L001_R1_001.fastq.gz
      [ID]_[LANE]_[READ]_[MULTIPLEX_ID].fastq.gz
      ID:CHX-022-T_S44
EOF


mkdir -p ${PDIR}/telomerecat
cd ${PDIR}/telomerecat/
mkdir input  output  scripts  sh_scripts

echo "${bold}WadingPool:${normal} Populating the input directory for Telomerecat directory...";
ln -s ${PDIR}/data/bam/*bam input/
ln -s ${PDIR}/data/bam/*bai input/
rm input/\*bam input/\*bai

cd sh_scripts
ln -s ${IDLIST} .
ln -s ${GIT}/preprocessing/telomerecat/generateScripts.sh .

echo "${bold}WadingPool:${normal} Running telomerecat...";
sh generateScripts.sh ${PDIR}
sh queueJobs.sh
