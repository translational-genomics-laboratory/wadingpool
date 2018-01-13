#!/bin/bash
#
#$ -cwd
#$ -S /bin/bash

module load vcftools/0.1.12
module load tabix/0.2.6
module load R/3.4.0


#### MERGE ALL VCFs TOGETHER
PDIR=$1
PDIR=${PDIR}'/variant_calling/mutect_common/output'
cd ${PDIR}

for eachDir in "filt_all filt_het"; do
#for eachDir in "filt_het"; do 
  echo ${eachDir}
  cd ${PDIR}/${eachDir}
  rm merge.vcf
  
  for i in $(ls . | grep "vcf$"); do
    bgzip $i
    tabix $i.gz
  done

  vcf-merge $(ls . | grep "vcf.gz$") > merge.vcf
done

