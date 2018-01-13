#!/bin/bash
#
#$ -cwd
#$ -S /bin/bash

module load vcftools/0.1.12
module load tabix/0.2.6
module load R/3.4.0

PDIR=$1
PDIR=${PDIR}'/variant_calling/mutect_common/output'
vcfFile=$2
vcfFile=${vcfFile}'.vcf'
#vcfFile='NET-2-011_02_T_DNA.processed.vcf'

cd ${PDIR}
mkdir filt_het filt_all


echo "Generating header file for ${vcfFile}"
grep "^#" chr1.${vcfFile} > ${vcfFile}.HEADER
HEADER=$(readlink -f ${vcfFile}.HEADER)

echo "Establishing the blank 'none' column..."
colidx=$(grep "#CHROM" chr1.${vcfFile} | \
  column -t |\
  perl -ne 'my @col = split(/\s+/, $_); if($col[9] =~ 'none') {print "10"} else {print "9"}')

cat $HEADER > filt_het/${vcfFile}
cat $HEADER > filt_all/${vcfFile}

for chr in $(seq 1 22); do
  echo "Adding on chromosome chr${chr}..."
  vcfid="chr${chr}.${vcfFile}"

  cat ${vcfid} | perl -ne 'my @col=split /\t/, $_; my @vcf=split /:/, $col['${colidx}']; if($vcf[3] > 1){ print $_ };' > filt_all/${vcfid}
  cat filt_all/${vcfid} >> filt_all/${vcfFile}

  cat filt_all/${vcfid} | perl -ne 'my @col=split /\t/, $_;
                           my @vcf=split /:/, $col['${colidx}'];
                           my @alleles=split /,/, $vcf[1];
                           if($alleles[0] != 0 && $alleles[1] != 0){ print $_};' >> filt_het/${vcfFile}
  rm filt_all/${vcfid}
done
