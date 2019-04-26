readCounter=$READCOPY
ichorScript=$ICHORPATH
ichorExtPath=/.mounts/labs/TGL/gsi/databases/ichorCNA_ext

mkdir -p $PDIR/ichor_cna/output/wig/

for ID in $(cat $IDLIST); do
  # link all bam files to tmp directory
  tmpDir=$PDIR/tmp_${ID}; mkdir -p $tmpDir
  ln -s $BAMDIR/$ID.bam.bam.cocleaned.bam $tmpDir/$ID.cocleaned.bam
  ln -s $BAMDIR/$ID.bam.bam.cocleaned.bai $tmpDir/$ID.cocleaned.bam.bai
  echo "module load R-gsi/3.5.1;
  $readCounter --window 1000000 --quality 20 \\
  --chromosome "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY" \\
  $tmpDir/$ID.cocleaned.bam | sed "s/chrom=chr/chrom=/" >  $PDIR/ichor_cna/output/wig/$ID.wig" > $ID.ichorCNA.sh
  ####
  echo "Rscript $ichorScript \\
  --id $ID \\
  --WIG $PDIR/ichor_cna/output/wig/$ID.wig \\
  --ploidy \"c(2,3)\" \\
  --normal \"c(0.2, 0.3, 0.4, 0.5,0.6,0.7,0.8,0.9)\" \\
  --maxCN 5 \\
  --gcWig $ichorExtPath/gc_hg19_1000kb.wig \\
  --mapWig $ichorExtPath/map_hg19_1000kb.wig\\
  --centromere $ichorExtPath/GRCh37.p13_centromere_UCSC-gapTable.txt \\
  --normalPanel $ichorExtPath/HD_ULP_PoN_1Mb_median_normAutosome_mapScoreFiltered_median.rds \\
  --includeHOMD False \\
  --chrs \"c(1:22, 'X')\" \\
  --chrTrain \"c(1:22)\" \\
  --estimateNormal True \\
  --estimatePloidy True \\
  --estimateScPrevalence True \\
  --scStates \"c(1,3)\" \\
  --txnE 0.9999 \\
  --txnStrength 10000 \\
  --outDir $PDIR/ichor_cna/output/ " >> $ID.ichorCNA.sh

  chmod +x $ID.ichorCNA.sh
  source $ID.ichorCNA.sh
done
