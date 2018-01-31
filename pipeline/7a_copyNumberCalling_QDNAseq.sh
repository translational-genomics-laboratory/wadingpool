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
ln -s ${GIT}/analysis/swgs/bin/QDNAseq_pipeline.R .

echo "${bold}WadingPool:${normal} Running QDNAseq copy-number caller...";
cat << EOF > runQDNAseq.sh
module load R/3.4.0

Rscript QDNAseq_pipeline.R \\
--pdir $PDIR \\
--outdir $RUNID \\
--runmode 'bin' \\
--binsize 50 \\
--regex $REGEX
EOF
qsub runQDNAseq.sh

