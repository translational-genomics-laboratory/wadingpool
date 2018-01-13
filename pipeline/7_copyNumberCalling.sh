## symlink FASTQ directory and align
bold=$(tput bold)
normal=$(tput sgr0)
cat << EOF
${bold}WadingPool:${normal} Variant calling 
    Helper script to process the aligned bam files and force variant calling using MuTect on set dbsnp_common sites
    Assumes that the bam file has been subsetted into chr-specific bams
      .
      ├── chr1
      │   ├── chr1.12_S4.bam
      │   └── chr1.12_S4.bam.bai
      ├── chr10
      │   ├── chr10.12_S4.bam
      │   └── chr10.12_S4.bam.bai
      ...
      └── chr9
          ├── chr9.12_S4.bam
          └── chr9.12_S4.bam.bai
EOF


mkdir -p ${PDIR}/qdnaseq
cd ${PDIR}/qdnaseq
mkdir bin_doc  input  output  sh_scripts

cp -R ${PDIR}/data/cocleaned_bam/symlinks/* input/

cd sh_scripts
ln -s ${IDLIST} .
ln -s ${GIT}/preprocessing/variant_calling/generateScripts.sh .

echo "${bold}WadingPool:${normal} Running telomerecat...";
sh generateScripts.sh ${PDIR} ${DBSNPDIR} ${DBSNPID}
sh queueJobs.sh

