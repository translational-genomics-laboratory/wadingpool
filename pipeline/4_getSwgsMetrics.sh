mkdir -p ${PDIR}/qc/picard
cd ${PDIR}/qc/picard
mkdir input  output_iSize  output_summary  output_wgs  scripts  sh_scripts

echo "${bold}WadingPool:${normal} Populating the input directory for QC directory...";
ln -s ${PDIR}/data/bam/*bam input/
ln -s ${PDIR}/data/bam/*bai input/
rm input/\*bam input/\*bai

cd sh_scripts
ln -s ${IDLIST} .
ln -s ${GIT}/preprocessing/qc_picard/generateScripts.sh .



echo "${bold}WadingPool:${normal} Subsetting each bam file by chromosome...";
sh generateScripts.sh ${PDIR}
sh queueJobs.sh