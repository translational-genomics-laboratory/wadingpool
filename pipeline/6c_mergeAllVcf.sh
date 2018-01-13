## symlink FASTQ directory and align
bold=$(tput bold)
normal=$(tput sgr0)


mkdir -p ${PDIR}/variant_calling/mutect_common/scripts
cd ${PDIR}/variant_calling/mutect_common/scripts

ln -s ${IDLIST} .
ln -s ${GIT}/preprocessing/variant_calling/mergeVcfs.sh .

echo "${bold}WadingPool:${normal} Merging all VCFs together into a single VCF...";
qsub mergeVcfs.sh ${PDIR}

