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
  PDIR='/mnt/work1/users/pughlab/projects/CHX/swgs/ichor_cna'
  readCounter='/mnt/work1/users/home2/quever/git/hmmcopy_utils/bin/readCounter'
  ichorPath='/mnt/work1/users/home2/quever/git/ichorCNA'
else
  echo "${bold}WadingPool:${normal} Using PDIR, dbsnp Directory and ID specified by configuration file"
  PDIR=$1
  readCounter=$2
  ichorPath=$3
fi
ichorDIR=${PDIR}'/ichor_cna/input'

rm queueJobs.sh
for ID in $(cat id_list.txt); do
  cat << EOF > ichorCNA.$ID.sh
#!/bin/bash
#
#$ -cwd

module load R/3.4.0

if [ ! -f $PDIR/input/$ID.cocleaned.bam.bai ]; then
  mv $PDIR/input/$ID.cocleaned.bai $PDIR/input/$ID.cocleaned.bam.bai
fi
$readCounter --window 1000000 --quality 20 \
--chromosome "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY" \
$PDIR/input/$ID.cocleaned.bam | sed "s/chrom=chr/chrom=/" >  $PDIR/output/wig/$ID.wig

Rscript $ichorPath/scripts/runIchorCNA.R \\
  --id $ID \\
  --WIG $PDIR/output/wig/$ID.wig \\
  --ploidy "c(2,3)" \\
  --normal "c(0.2, 0.3, 0.4, 0.5,0.6,0.7,0.8,0.9)" \\
  --maxCN 5 \\
  --gcWig $ichorPath"/inst/extdata/gc_hg19_1000kb.wig" \\
  --mapWig $ichorPath"/inst/extdata/map_hg19_1000kb.wig" \\
  --centromere $ichorPath"/inst/extdata/GRCh37.p13_centromere_UCSC-gapTable.txt" \\
  --normalPanel $ichorPath"/inst/extdata/HD_ULP_PoN_1Mb_median_normAutosome_mapScoreFiltered_median.rds" \\
  --includeHOMD False --chrs "c(1:22, \\"X\\")" --chrTrain "c(1:22)" \\
  --estimateNormal True --estimatePloidy True --estimateScPrevalence True \\
  --scStates "c(1,3)" --txnE 0.9999 --txnStrength 10000 --outDir $PDIR/output/
  

EOF

  
  echo "qsub ichorCNA."$ID".sh" >> queueJobs.sh
done