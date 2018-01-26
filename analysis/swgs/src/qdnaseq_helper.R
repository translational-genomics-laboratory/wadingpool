# Run the default QDNAseq pipeline for a given bin-size and set of bam files
getDefaultCn <- function(bins, bam.files, bam.dir, run.em=TRUE){
  filt.idx <- which(pData(bins)$end < pData(bins)$start)
  if(length(filt.idx) > 0L){
    bins <- bins[-filt.idx,]
  }
  # Bins reads into each of the Xkb size bins
  readCounts <- binReadCounts(bins, 
                              bamfiles=bam.files,
                              path=bam.dir, 
                              isPaired=TRUE, isProperPair=TRUE,
                              isUnmappedQuery=FALSE, hasUnmappedMate=FALSE,
                              isNotPassingQualityControls=FALSE, isDuplicate=FALSE)
  readCountsFiltered <- applyFilters(readCounts,
                                     residual=TRUE, blacklist=TRUE) # 1000genomes blacklist
  
  # Fits a surface Loess to (mappability x GC) matrix to correct each bin by
  readCountsFiltered <- estimateCorrection(readCountsFiltered)
  # Corrects the read counts by dividing by the fitted Loess (mappability x GC) matrix
  copyNumbers <- correctBins(readCountsFiltered, method="ratio")  #ratio (default) divides counts with fit
  copyNumbersNormalized <- normalizeBins(copyNumbers, force=TRUE) # median centers
  copyNumbersSmooth <- smoothOutlierBins(copyNumbersNormalized)
  
  #  Segments using CBS following an Anscombe transformation
  copyNumbersSegmented <- segmentBins(copyNumbersSmooth, transformFun="sqrt")
  # "Normalizes" by log2 variance stabilizing, then centering at highest density peak, then reverse the transformation
  copyNumbersSegmented <- normalizeSegmentedBins(copyNumbersSegmented)
  
  if(run.em){
    # Applies the 6-mer gaussian deconvolution to assign copy-state used in CGHcall
    copyNumbersCalled <- callBins(copyNumbersSegmented)  
  } else {
    copyNumbersCalled <- NULL
  }
  list("copyNumbersCalled" = copyNumbersCalled,
       "copyNumbersSegmented" = copyNumbersSegmented,
       "copyNumbersSmooth" = copyNumbersSmooth,
       "copyNumbersNormalized" = copyNumbersNormalized,
       "copyNumbers" = copyNumbers,
       "readCountsFiltered" = readCountsFiltered,
       "readCounts"=readCounts)
}

# For a set of bins that overlap, intersect get the mapping indices
overlapQdnaBins <- function(x, bins, cn.type="copynumber"){
  print(head(pData(featureData(x)), n=1))
  isec.idx <- interesctAdf(bins, x)
  if(cn.type == 'copynumber'){
    cnvals  <- QDNAseq:::copynumber(x)[isec.idx,]
  } else if(cn.type =='segmented'){
    cnvals  <- QDNAseq:::segmented(x)[isec.idx,]
  }
  cnvals 
}

# Using the intersect/mapping indices, get the estimates of mean + sd
qdnaMeanSd <- function(cnValues, each.id){
  col.idx <- grep(each.id, colnames(cnValues))
  summ.cn <- t(apply(cnValues[,col.idx], 1, function(x) {
    m <- mean(x, na.rm = TRUE)
    s <- sd(x, na.rm = TRUE)
    c("mean"=m, "sd"=s)
  }))
  summ.cn[is.nan(summ.cn[,'mean']), 'mean'] <- NA
  summ.cn
}

# Function to generate the CBS segments
getCbsSeg <- function(chr, val, pos, name, aval=0.01, nperm=1000, weights=NULL){
  chr.id <- unlist(regmatches(chr, regexec(pattern = "[0-9X]+", chr)))
  if(!chr.id %in% 'X'){
    chr.num <- as.numeric(chr.id)
  } else {
    chr.num <- chr.id
  }
  
  CNA.object <- CNA(as.numeric(val),
                    rep(chr.num, length(pos)),
                    pos, data.type="logratio", sampleid=name)
  smoothed.CNA.object <- smooth.CNA(CNA.object)
  if(!is.null(weights)){
    if(any(is.na(weights))) weights[is.na(weights)] <- 0.0001
    if(any(weights == 0)) weights[weights == 0] <- 0.0001
    segment.smoothed.CNA.object <- segment(CNA.object, verbose=1, 
                                           alpha=aval, nperm=nperm,
                                           weights=weights)
  } else {
    segment.smoothed.CNA.object <- segment(CNA.object, verbose=1, 
                                           alpha=aval, nperm=nperm)
  }
  return(segment.smoothed.CNA.object)
}

#  Generate the QDNAseq results
plotResults <- function(pdir, type, ...){
  cat("Plotting results:
      \t'all' = Generate all plots
      \t'raw' = Generate raw and QC plots
      \t'smooth' = Generate mappability/GC corrected raw plots
      \t'segmented' = Generate the CBS segmented plots
      \t'called' = Generate the 6mer EM-deconstructed CN Calls")

  if(type == 'raw' || type == 'all'){
    cat("Plotting Raw and QC plots...")
    pdf(file.path(pdir, "qdnaseq", "output", opt$outdir, "readCounts.RawAndQC.pdf"))
    isobarPlot(readCountsFiltered, what="fit")
    noisePlot(readCountsFiltered)
    plot(readCounts, logTransform=FALSE, ylim=c(quantile(readCounts@assayData$counts, 0.01),
                                                quantile(readCounts@assayData$counts, 0.99)))
    try(expr=highlightFilters(readCounts, logTransform=FALSE,
                              residual=TRUE, blacklist=TRUE), 
        silent = TRUE)
    dev.off()
  }
  
  if(type == 'smooth' || type == 'all'){
    cat("Plotting Mappability-GC smoothed copy-ratio plots...")
    pdf(file.path(pdir, "qdnaseq", "output", opt$outdir, "readCounts.smooth.pdf"))
    plot(copyNumbersSmooth, ylim=c(0,15))
    dev.off()
  }
  
  if(type == 'segmented' || type == 'all'){
    cat("Plotting CBS segmented plots...")
    pdf(file.path(pdir, "qdnaseq", "output", opt$outdir, "readCounts.segmented.pdf"))
    plot(copyNumbersSegmented)
    dev.off()
  }
  
  if(type == 'called' || type == 'all'){
    cat("Plotting Copy-states plots...")
    pdf(file.path(pdir, "qdnaseq", "output", opt$outdir, "readCounts.called.pdf"))
    plot(copyNumbersCalled)
    dev.off()
  }
  
  if(opt$devel){
    pdf(file.path(pdir, "qdnaseq", "output", opt$outdir, "readCounts.ascn.pdf"))
    apply(copyNumbersASCN, 2, function(x) plot(x, ylim=c(0,6), pch=15))
    dev.off()
  }
  
  if(type == 'segmented' ||type == 'all'){
    if(type == 'segmented'){
      segdat <- formatSeg(copyNumbersSegmented)
    } else {
      segdat <- formatSeg(copyNumbersCalled, method='calls')
      save(copyNumbersSegmented, copyNumbersCalled, 
           bins, file=file.path(pdir, "qdnaseq", "output", opt$outdir, "readCounts.Rdata"))
      print(file.path(pdir, "qdnaseq", "output", opt$outdir, "readCounts.Rdata"))
    }
    
    write.table(do.call("rbind", segdat),
                file=file.path(pdir, "qdnaseq", "output", opt$outdir, "swgs_segments.seg"),
                col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
    print(file.path(pdir, "qdnaseq", "output", opt$outdir, "swgs_segments.seg"))
    
  }
}

argChecker <- function(opt, opt_parser){
  if (is.null(opt$pdir)){
    print_help(opt_parser)
    stop("Parent directory must be supplied", call.=FALSE)
  } else if(is.null(opt$regex)){
    print("WARNING: No regex was given so all bam files in the input directory will be loaded")
  } else if(opt$mode == 'stagger' && opt$binsize != 1000){
    print("WARNING: In stagger mode, all bin-sizes are default to 1000kb bins with 5kb offsets.  Input binsize will be ignored.")
  }
}

