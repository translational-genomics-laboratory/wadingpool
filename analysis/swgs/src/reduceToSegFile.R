formatSeg <- function(cnseg, method="segmented"){
  fdat <- Biobase::fData(cnseg)
  seginf <- QDNAseq:::segmented(cnseg)
  if (method=="calls") { callinf <- QDNAseq:::calls(cnseg) }
  
  appendSegToFdata <- function(x, ...){
    idval <- gsub(".cocleaned$", "", rep(x, nrow(seginf)))
    chrpos <- c("chromosome", "start", "end")
    comb_segdat <- cbind(idval, 
                         fdat[,chrpos],
                         seginf[,grep(x, colnames(seginf)), drop=FALSE],
                         fdat[,which(!colnames(fdat) %in% chrpos)])
    colnames(comb_segdat)[1:5] <- c("ID", "chrom", "loc.start", "loc.end", "seg.mean")
    
    if(method=='calls') {
      comb_segdat <- cbind(comb_segdat, callinf[,grep(x, colnames(callinf)), drop=FALSE] + 2)
      colnames(comb_segdat)[ncol(comb_segdat)] <- 'cn'
    }
    comb_segdat
  }
  comb_seg <- lapply(colnames(seginf), appendSegToFdata)
  names(comb_seg) <- gsub(".cocleaned$", "", colnames(seginf))
  
  reduceSeg <- function(x, ...){
    red_chr <- lapply(split(x, f=x$chrom), function(each_chr){
      tryCatch({
        each_chr <- each_chr[-which(is.na(each_chr$seg.mean)),]
        rle_idx <- getRleIdx(each_chr, "seg.mean")
        
        cn_seg <- data.frame("ID"=unique(each_chr$ID),
                             "chrom"=unique(each_chr$chrom),
                             "loc.start"=each_chr$loc.start[rle_idx$start.idx],
                             "loc.end"=each_chr$loc.end[rle_idx$end.idx],
                             "num.mark"=rle_idx$lengths,
                             "seg.mean"=round(as.numeric(as.character(rle_idx$values)),4))
        if(method == 'calls'){
          cnval <- sapply(seq(length(rle_idx$start.idx)), function(ridx){
            mean(each_chr[rle_idx$start.idx[ridx]:rle_idx$end.idx[ridx], 'cn'], na.rm=TRUE)
          })
          
          cn_seg <- cbind(cn_seg, "CN"=cnval)
        }
        cn_seg
      }, error=function(e){NA})
    })
    red_chr <- as.data.frame(do.call("rbind", red_chr[as.character(c(1:22, "X", "Y"))]))
    red_chr <- red_chr[-which(is.na(red_chr$ID)),]
    red_chr$seg.mean <- log2(red_chr$seg.mean)
    red_chr
  }
  red_seg <-  lapply(comb_seg, reduceSeg)

  red_seg
}