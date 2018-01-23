### try http:// if https:// URLs are not supported
# source("http://bioconductor.org/biocLite.R")
# biocLite("QDNAseq")
# 
#  Rscript QDNAseq_pipeline.R /path/to/parent_directory outdir_id run_mode bin_size BAM_004_076 
#     args[1] = Absolute path to the parent directory;  Assumes nested directories are 'qdnaseq' and 'input'
#     args[2] = Name of the out-directory (i.e. 'CHX_batch1')
#     args[3] = Mode to run in: 'stagger' or 'bin'
#     args[4] = Bin size to use:  automatically defaults to 1000kb if using 'stagger' mode
#     args[5] = [OPTIONAL]: Regex to select for specific BAM files

library(QDNAseq)
library(GenomicRanges)
library(Biobase)
library(DNAcopy)
source("~/git/wadingpool/analysis/swgs/src/misc.R")
source("~/git/wadingpool/analysis/swgs/src/physicalCov.R")
source("~/git/wadingpool/analysis/swgs/src/reduceToSegFile.R")
source("~/git/wadingpool/analysis/swgs/src/qdnaseq_helper.R")

#### Initial setup
########################################
printUsage()

args = commandArgs(trailingOnly=TRUE)
devel <- FALSE  # Enable development mode (Default = FALSE)
if (length(args)==0) {
  pdir <- '/mnt/work1/users/pughlab/projects/Tsao_lung_organoid/'
  pdir <- '/mnt/work1/users/pughlab/projects/NET-SEQ/shallow_wgs'
  pdir <- '/mnt/work1/users/pughlab/projects/shallowWGS_benchmark/myeloma_cfdna/signy-myelstone'
  outdir <- "004-076"
  run.mode <- 'stagger'
  bin.size <- 1000
  bam.grep <- "004-076" # NULL
} else {
  pdir <- args[1]
  outdir <- args[2]
  run.mode <- args[3]
  if(run.mode == 'stagger') bin.size <- 1000 else bin.size <- as.integer(as.character(args[4]))
  if(length(args)==5L) bam.grep <- args[5] else bam.grep <- NULL
}
offset.dir <- '~/git/wadingpool/analysis/swgs/ref/stagger'
bam.dir <- file.path(pdir, "qdnaseq", "input")
setwd(bam.dir)


#### Copy-number calling Pipeline set up by QDNAseq 
########################################
# Adjust what bam files are included in the analysis
bam.files <- list.files(bam.dir, pattern="bam$")
if(!is.null(bam.grep)) bam.files <- bam.files[grep(bam.grep, bam.files)]

if(bin.size == 1000){
  load(file.path(offset.dir, "1000kb.0.Rda"))  # Default with no offset/staggering
} else {
  bins <- getBinAnnotations(binSize=bin.size)  
}

# Set up the default copy-number estimates 
cnList <- getDefaultCn(bins, bam.files, bam.dir, run.em=TRUE)

if(run.mode == 'stagger'){
  copyNumbersSegmented <- cnList[['copyNumbersSegmented']]
  # Set up a reference 5kb copy-number estimates to refine
  bins.5kb <- getBinAnnotations(binSize=5)  #per Kb5
  bins.5kb <- getDefaultCn(bins.5kb, bam.files, bam.dir, run.em=FALSE)[['copyNumbersSegmented']]
  
  # For each staggered bin, get the copy-number estimates
  window.copyNumbersSegmented <- lapply(list.files(offset.dir), function(binx){
    load(file.path(offset.dir, binx))
    getDefaultCn(bins, bam.files, bam.dir, run.em=FALSE)[['copyNumbersSegmented']]
  })
  names(window.copyNumbersSegmented) <- list.files(offset.dir)
  
  
  # Intersect and take a mean+sd estimate around the copy-number ratios and the segment values
  cnValues <- lapply(window.copyNumbersSegmented, overlapQdnaBins, 
                     bins=bins.5kb, cn.type='copynumber')
  cnValues <- do.call("cbind", cnValues)
  segValues <- lapply(window.copyNumbersSegmented, overlapQdnaBins, 
                      bins=bins.5kb, cn.type='segmented')
  segValues <- do.call("cbind", segValues)
  
  # For copy-number or segmented, estimate the mean and sd for each list
  summ.list <- list()
  cn.type <- 'copynumber' #segmented, copynumber
  for(each.id in colnames(QDNAseq:::copynumber(bins.5kb))){
    if(cn.type == 'segmented'){
      summ.cn <- qdnaMeanSd(segValues, each.id)
      QDNAseq:::segmented(bins.5kb)[,each.id] <- summ.cn[,1]
    } else if(cn.type == 'copynumber'){
      summ.cn <- qdnaMeanSd(cnValues, each.id)
      QDNAseq:::copynumber(bins.5kb)[,each.id] <- summ.cn[,1]
    }
    
    summ.list[[each.id]] <- summ.cn
  }
  
} else {
  copyNumbersCalled <- cnList[["copyNumbersCalled"]]
  copyNumbersSegmented <- cnList[["copyNumbersSegmented"]]
  copyNumbersSmooth <- cnList[["copyNumbersSmooth"]]
  copyNumbersNormalized <- cnList[["copyNumbersNormalized"]]
  copyNumbers <- cnList[["copyNumbers"]]
  readCountsFiltered <- cnList[["readCountsFiltered"]]
  readCounts <- cnList[["readCounts"]]
}

if(devel){
  # Plot the means with a standard deviation, resegment using a SD-weighted CBS
  dir.create(file.path(pdir, "qdnaseq", "output"), recursive = TRUE)
  pdf(file.path(pdir, "qdnaseq", "output", "trial.pdf"), width = 14)
  for(each.file in names(summ.list)){
    sf <- summ.list[[each.file]]
    weights <- (1- (sf[,'sd'] / quantile(sf[,'sd'], 0.95, na.rm = TRUE)))
    weights[weights< 0] <- 0
    na.vals <- which(is.na(sf[,'mean']))
    
    sf.cbs <- getCbsSeg(chr='1', val=(sf[-na.vals,'mean']),
                        pos = c(1:nrow(sf))[-na.vals], name="test",
                        weights=weights[-na.vals],
                        aval=0.05, nperm=1000)
    sf.seg <- sf.cbs$output
    
    
    plot(c(1:nrow(sf)), sf[,'mean'], type='n', main=each.file)
    polygon(c(1:nrow(sf), nrow(sf):1), c((sf[,'mean'] + sf[,'sd']), 
                                         rev(sf[,'mean'] - sf[,'sd'])),
            col = "grey30", border = NA)
    for(i in seq(1:dim(sf.seg)[1])) {
      #guideline at CNV median
      segments(sf.seg[i,'loc.start'],
               as.numeric(sf.seg[i,'seg.mean']),
               sf.seg[i,'loc.end'],
               as.numeric(sf.seg[i,'seg.mean']),
               col="orange",
               lwd=2)       
    }  
    #points(c(1:nrow(sf)), sf[,'mean'], col="black", pch=20)
  }
  dev.off()
}




#### QDNAseq Plotting (temporary until heterozygostiy is included)
########################################
dir.create(file.path(pdir, "qdnaseq", "output", outdir), recursive = TRUE)
cat("Saving image...\n")
save.image(file.path(pdir, "qdnaseq", "output", outdir, "CHXbatch3.Rdata"))
cat("Plotting...\n")
plotResults(pdir, type="all")   #raw, smooth, segmented, called



####
if(devel){
  fitpro <- getEMestimates(copyNumbersSegmented, break.num=500)
  seg_ascn <- copyNumbersSegmented@assayData$segmented
  for(each_sample in names(fitpro)){
    seg_ascn[,each_sample] <- as.integer(as.character(fitpro[[each_sample]]$prob.cn.df))
  }
  
  as.integer(fitpro[[1]]$prob.cn.df)
  fitpro <- getEMestimates(wgsCov.df, break.num=800)
  
  
  dir.create(file.path(pdir, "bin_doc", "output"), recursive = TRUE)                   
  pdf(file.path(pdir, "bin_doc", "output", "readCounts.raw.pdf"))
  isobarPlot(readCountsFiltered, what="fit")
  #noisePlot(readCountsFiltered)
  
  plot(readCounts, logTransform=FALSE, ylim=c(quantile(readCounts@assayData$counts, 0.01), 
                                              quantile(readCounts@assayData$counts, 0.99)))
  highlightFilters(readCounts, logTransform=FALSE,
                   residual=TRUE, blacklist=TRUE)
  plot(copyNumbersSmooth, ylim=c(0,15))
  plot(copyNumbersSegmented)
  #plot(copyNumbersCalled)
  
  dev.off()
  
  
}
