mkdir -p ${PDIR}/qc/picard
cd ${PDIR}/qc/picard
mkdir input  output_iSize  output_summary  output_wgs  scripts  sh_scripts

echo "${bold}WadingPool:${normal} Populating the input directory for QC directory...";
cp -R ${PDIR}/data/cocleaned_bam/symlinks/* input/

cd sh_scripts
rm id_list.txt generateScripts.sh
ln -s ${IDLIST} .
ln -s ${GIT}/preprocessing/qc_picard/generateScripts.sh .



echo "${bold}WadingPool:${normal} Subsetting each bam file by chromosome...";
sh generateScripts.sh ${PDIR}
sh queueJobs.sh