## symlink FASTQ directory and align
bold=$(tput bold)
normal=$(tput sgr0)


mkdir -p ${PDIR}/variant_calling/mutect_common/scripts
cd ${PDIR}/variant_calling/mutect_common/scripts

ln -s ${IDLIST} .
ln -s ${GIT}/preprocessing/variant_calling/getHetSnps.sh .

echo "${bold}WadingPool:${normal} filtering vcf for only heterozygous and covered variants...";
cat << EOF > queueJobs.sh
for id in $(cat id_list.txt); do
  qsub getHetSnps.sh ${PDIR} ${id}
done
EOF
sh queueJobs.sh

