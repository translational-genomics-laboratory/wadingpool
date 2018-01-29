mkdir -p ${PDIR}/data/sh_scripts/chr_subset
cd ${PDIR}/data/sh_scripts/chr_subset

rm id_list.txt generateScripts.sh 
ln -s ${IDLIST} .
ln -s ${GIT}/preprocessing/alignment/chr_subset/generateScripts.sh .


echo "${bold}WadingPool:${normal} Subsetting each bam file by chromosome...";
sh generateScripts.sh ${PDIR}
sh queueJobs.sh