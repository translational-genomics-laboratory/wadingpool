suppressPackageStartupMessages(library("QDNAseq"))
suppressPackageStartupMessages(library("GenomicRanges"))
suppressPackageStartupMessages(library("Biobase"))
suppressPackageStartupMessages(library("DNAcopy"))
suppressPackageStartupMessages(library("optparse"))
source("~/git/wadingpool/analysis/swgs/src/misc.R")
source("~/git/wadingpool/analysis/swgs/src/physicalCov.R")
source("~/git/wadingpool/analysis/swgs/src/reduceToSegFile.R")
source("~/git/wadingpool/analysis/swgs/src/qdnaseq_helper.R")
offset.dir <- '~/git/wadingpool/analysis/swgs/ref/stagger'

#### Initial setup
########################################
args = commandArgs(trailingOnly=TRUE)

option_list = list(
  make_option(c("-p", "--pdir"), type="character", default=NULL, 
              help="Absolute path to the Parent DIRectory which contains the qdnaseq directory", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default="outdir", 
              help="Name of the out-directory that you want to save your results to [default= %default]", metavar="character"),
  make_option(c("-m", "--mode"), type="character", default="bin", 
              help="Copy-number calling for 'bin' or 'stagger' mode [default= %default]", metavar="character"),
  make_option(c("-b", "--binsize"), type="integer", default=1000, 
              help="Size of the genomic bins [default= %default]", metavar="integer"),
  make_option(c("-r", "--regex"), type="character", default=NULL, 
              help="[OPTIONAL] Regex to grep certain bam files for selective processing [default= %default]", metavar="character"),
  make_option(c("-d", "--devel"), type="logical", default=FALSE, 
              help="Developer mode [default= %default]", metavar="logical")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
argChecker(opt, opt_parser)
#print_help(opt_parser)


if(0){
  # Used for running in interactive mode
  opt$pdir <- '/mnt/work1/users/pughlab/projects/Tsao_lung_organoid/'
  opt$pdir <- '/mnt/work1/users/pughlab/projects/NET-SEQ/shallow_wgs'
  opt$pdir <- '/mnt/work1/users/pughlab/projects/shallowWGS_benchmark/myeloma_cfdna/signy-myelstone'
  opt$outdir <- "004-076"
  opt$runmode <- 'bin'
  opt$binsize <- 1000
  opt$regex <- "004-076" # NULL
} 
bam.dir <- file.path(opt$pdir, "qdnaseq", "input")
setwd(bam.dir)



#### Copy-number calling Pipeline set up by QDNAseq 
########################################
# Adjust what bam files are included in the analysis
bam.files <- list.files(bam.dir, pattern="bam$")
if(!is.null(opt$regex)) bam.files <- bam.files[grep(opt$regex, bam.files)]

if(opt$binsize == 1000){
  load(file.path(offset.dir, "1000kb.0.Rda"))  # Default with no offset/staggering
} else {
  bins <- getBinAnnotations(binSize=opt$binsize)  
}

# Set up the default copy-number estimates 
cnList <- getDefaultCn(bins, bam.files, bam.dir, run.em=TRUE)

if(opt$runmode == 'stagger'){
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

if(opt$devel){
  # Plot the means with a standard deviation, resegment using a SD-weighted CBS
  dir.create(file.path(opt$pdir, "qdnaseq", "output"), recursive = TRUE)
  pdf(file.path(opt$pdir, "qdnaseq", "output", "trial.pdf"), width = 14)
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
dir.create(file.path(opt$pdir, "qdnaseq", "output", opt$outdir), recursive = TRUE)
cat("Saving image...\n")
save.image(file.path(opt$pdir, "qdnaseq", "output", opt$outdir, "CHXbatch3.Rdata"))
cat("Plotting...\n")
plotResults(opt$pdir, type="all")   #raw, smooth, segmented, called



####
if(opt$devel){
  fitpro <- getEMestimates(copyNumbersSegmented, break.num=500)
  seg_ascn <- copyNumbersSegmented@assayData$segmented
  for(each_sample in names(fitpro)){
    seg_ascn[,each_sample] <- as.integer(as.character(fitpro[[each_sample]]$prob.cn.df))
  }
  
  as.integer(fitpro[[1]]$prob.cn.df)
  fitpro <- getEMestimates(wgsCov.df, break.num=800)
  
  
  dir.create(file.path(opt$pdir, "bin_doc", "output"), recursive = TRUE)                   
  pdf(file.path(opt$pdir, "bin_doc", "output", "readCounts.raw.pdf"))
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
