#!/bin/bash
#
#$ -cwd
#$ -S /bin/bash
#
############################################################
bold=$(tput bold)
normal=$(tput sgr0)

## Specify PDIR and LIBID as passed in or hardcoded
if [ $# -eq 0 ]; then
  echo "${bold}WadingPool:${normal} Using hardcoded PDIR, dbsnp Directory and ID"
  PDIR='/mnt/work1/users/pughlab/projects/CHX/swgs'
  dbsnpDir='/mnt/work1/users/pughlab/references/dbsnp/dbsnp_b146_GRCh37p13/common/bed'
  dbsnpId='common_all_20151104.bed'
else
  echo "${bold}WadingPool:${normal} Using PDIR, dbsnp Directory and ID specified by configuration file"
  PDIR=$1
  dbsnpDir=$2
  dbsnpId=$3
fi
PDIR=${PDIR}'/variant_calling/mutect_common'

rm queueJobs.sh
for chr in $(seq 1 22); do
  for id in $(cat id_list.txt); do
    cat << EOF > mutect.${chr}.${id}.sh
#!/bin/bash
#
#$ -cwd

TUMOR="${id}"
dbsnpDir=${dbsnpDir}
dbsnpId=${dbsnpId}
PDIR=${PDIR}

module load igenome-human/hg19
module load  mutect/1.1.5
module load gatk/3.5
module unload java/6
module load java/7

java -Xmx13g -jar \$mutect_dir/muTect.jar \\
--analysis_type  MuTect \\
--min_qscore 20 \\
-R \${REF} \\
-I:tumor \${PDIR}/input/chr_subset/chr${chr}/chr${chr}.\${TUMOR}.bam \\
-L \${dbsnpDir}/chr${chr}.\${dbsnpId} \\
--force_output \\
-vcf \${PDIR}/output/chr${chr}.\${TUMOR}.vcf

EOF
  
    echo "qsub mutect.${chr}.${id}.sh" >> queueJobs.sh
  done
done