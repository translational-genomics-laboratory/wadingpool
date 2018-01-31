## symlink FASTQ directory and align
bold=$(tput bold)
normal=$(tput sgr0)
cd ${PDIR}/telomerecat/sh_scripts
rm id_list.txt runTelolen.sh
ln -s ${IDLIST} .
ln -s ${GIT}/preprocessing/telomerecat/runTelolen.sh .

echo "${bold}WadingPool:${normal} Running telomerecat...";
qsub runTelolen.sh ${PDIR}/telomerecat

