library(markovchain)
library(optparse)
library(Hmisc)

option_list=list(
  make_option(c("-i", "--id"), type='character', default=NULL,
              help='ID of the file you are passing in [id].pos,.tags.list,.ld.gz', metavar='character'),
  make_option(c("-d", "--dir"), type='character', default=NULL,
              help='Absolute path to directory that the input and output files are to be written to', metavar='character'),
  make_option(c("-m", "--mode"), type='character', default='model',
              help='Mode to run the script in.  "model" will calculate non-homogenous markov-models and binomial probabilities and save them. "sample" will sample from the MC and Binom models',
              metavar='character'),
  make_option(c("-r", "--ref"), type='character', default='0',
              help='State name of REF/Haplotype_1 allele',
              metavar='character'),
  make_option(c("-a", "--alt"), type='character', default='1',
              help='State name of ALT/Haplotype_2 allele',
              metavar='character'),
  make_option(c("-w", "--weight"), type='integer', default='100',
              help='Weight of importance for founder SNPs to sub-founder SNPs',
              metavar='character'),
  make_option(c("-G", "--geno"), type='character', default=NULL,
              help='OPTIONAL: If running all haploblocks, use this in pair with --genoid.  geno=name of the VCF',
              metavar='character'),
  make_option(c("-g", "--genoid"), type='character', default=NULL,
              help='OPTIONAL: If running all haploblocks, use this in pair with --geno.  genoid=ID associated with VCF',
              metavar='character')
)

opt_parser=OptionParser(usage = "Rscript %prog --input <INPUT FILE>
                        --output <OUTPUT FILE>
                        --path <OUTPUT DIRECTORY>", option_list=option_list)
opt=parse_args(opt_parser)

#########
### FUNCTIONS

## Modelling functions
genSmallTm <- function(R){
  R <- R/2 + 0.5
  matrix(c(R, 1-R, 1-R, R), ncol=2)
}

generateNHMC <- function(marker, states, prob){
  fMC <- lapply(seq_along(marker[-1]), function(i){
    new("markovchain", states = states,
        transitionMatrix = genSmallTm(prob[i]),
        name = paste0("state t", (i-1)))
  })
  mcCCRC <- new("markovchainList", markovchains = fMC,
                name = "Founder SNPs")
  mcCCRC
}

generateBinProb <- function(x, founders.idx, mat){
  subfounders.R <- mat[founders.idx[x], c((founders.idx[x]+1):(founders.idx[x+1]-1)), drop=FALSE]
  subfounders.R <- (subfounders.R / 2) + 0.5
  subfounders <-  colnames(subfounders.R)
  
  subfounders.R
}

filterPos <- function(pos.df, all.ids, fidx){
  ids.in.block <- all.ids[min(fidx):max(fidx)]
  #pos.df$V3 <- gsub(";.*", "", pos.df$V3)
  pos.df <- pos.df[which(pos.df$V3 %in% ids.in.block), ]
  pos.df
}

## Sampling functions
sampleFromBinomial <- function(x, n){
  r.match <- rbinom(n=(length(x) * n), size=1, prob=as.numeric(x))
  r.match <- matrix(r.match, ncol=length(x), byrow=TRUE)
  r.match
}

sampleFromNHMC <- function(model, init.prob, n, states){
  m.sim <- lapply(rbinom(n=n, size=1, prob=init.prob), function(init.state){
    init.state <- states[(init.state+1)]
    as.character(rmarkovchain(n = 1, object = model, t0 = init.state, include.t0 = TRUE)$values)
  })
  m.sim
}

combineSampling <- function(mc.s, states, bin.s, ids){
  all.sim<- lapply(mc.s, function(f.sim){
    f.with.subf <- lapply(seq_along(f.sim[-length(f.sim)]), function(each.f) {
      mismatch.f <- states[states != f.sim[each.f]]
      subf.hap <- matrix(mismatch.f, ncol=ncol(bin.s[[each.f]]), nrow=nrow(bin.s[[each.f]]))
      subf.hap[bin.s[[each.f]] == 1] <- f.sim[each.f]
      
      cbind(f.sim[each.f], subf.hap)
    })
    f.with.subf[[length(f.sim)]] <- matrix(rep(f.sim[length(f.sim)], opt$weight), ncol=1)
    f.with.subf <- do.call("cbind", f.with.subf)
  })
  all.sim <- do.call("rbind", all.sim)
  if(all.is.numeric(all.sim)) storage.mode(all.sim) <- "numeric"
  colnames(all.sim) <- ids
  
  all.sim
}

###### Main Function
LDMC <- function(each.row){
  err <- NULL
  
  ##/ Start of helper function
  if(RUNHELPER){
    opt$id <- blocks.det$range[each.row]
    founders.mat <- as.matrix(strsplit(blocks.det[each.row, 'SNPS'], split="\\|")[[1]])
    
    # Creating the founders table
    write.table(founders.mat, file=file.path(opt$dir, paste0(opt$id, "_founders.txt")),
                quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE)
    
    # Running the helper shell script for plink
    command <- paste("sh", file.path(base.dir, "helperLdMat.sh"), 
                     opt$geno, opt$genoid, opt$id, sep=" ")
    
    setwd(base.dir)
    system(command, intern = FALSE,
           ignore.stdout = TRUE, ignore.stderr = FALSE,
           wait = TRUE)
    setwd(opt$dir)
  }
  ##/ End of helper function
  
  
  ##  Setting up File IDs
  pos.id <- paste0(opt$id, ".pos")
  tags.id <- paste0(opt$id, ".tags.list")
  maf.id <- paste0(opt$id, ".maf.ld.gz")
  ld.mat <- paste0(opt$id, ".ld.gz")
  founders.ids <- paste0(opt$id, "_founders.txt")
  mcmodel.id <- paste0(opt$id, "_mcModel.rda")
  
  ## Modelling:
  if(grepl("model", opt$mode, ignore.case=TRUE)){
    ## Validate files exist:
    exists.flag <- !all(file.exists(paste0(opt$dir, 
                                           paste0("/", opt$id, 
                                                  c(".log", ".maf.log", ".nosex", 
                                                    ".tags.list", "_founders.txt", 
                                                    ".ld.gz", ".maf.ld.gz", ".pos", ".maf.nosex")))))
    if(exists.flag) {
      err <- each.row
      next(paste0("Could not find all files for haploblock", each.row))
    }
    
    ##  Reading in files
    pos.df <- read.table(pos.id, header=FALSE, stringsAsFactors=FALSE)
    founders <- read.table(founders.ids, header=FALSE, stringsAsFactors=FALSE)$V1
    tags.df <- read.table(tags.id, header=TRUE, stringsAsFactors=FALSE)
    maf.df <- read.table(maf.id, header=TRUE, stringsAsFactors=FALSE)
    ld.mat <- round(read.table(ld.mat, header=FALSE), 3)
    
    if(any(duplicated(tags.df$SNP))) {
      err <- each.row
      dup.idx <- which(duplicated(tags.df$SNP))
      warning(paste0("Duplicates found for Block: ", each.row, ", Marker ID: ", tags.df$SNP[dup.idx], "\n"))
      ld.mat <- ld.mat[-dup.idx, -dup.idx,drop=FALSE]
      tags.df <- tags.df[-dup.idx,,drop=FALSE]
      pos.df <- pos.df[-dup.idx, ,drop=FALSE]
    }
    #colnames(ld.mat) <- rownames(ld.mat) <- gsub(";.*", "", tags.df$SNP)
    colnames(ld.mat) <- rownames(ld.mat) <- tags.df$SNP
    
    if(!all(founders %in% tags.df$SNP)) {
      err <- each.row
      next(paste0("Founders didn't match the SNPs, likely failure in plink run: haploblock ", each.row))
    }
    
    
    # Calculate step-wise R^2 between founders
    #founders <- as.character(sapply(founders, function(x) gsub(";.*", "", x)))
    founders <- as.character(sapply(founders, function(x) x))
    founders.idx <- sapply(founders, function(x) which(x == colnames(ld.mat)))
    if(is.list(founders.idx)){
      rm.idx <- which(sapply(founders.idx, length) == 0)
      warning(paste0("Warning: could not find the founder SNP, ", paste(founders[rm.idx], collapse='/'), ", in the LD matrix. Removing..."))
      founders <- founders[rm.idx]
      founders.idx <- as.integer(founders.idx[-rm.idx])
    }
    founders.R2 <- diag(as.matrix(ld.mat[founders.idx[-length(founders.idx)], founders.idx[-1]]))
    if(any(is.nan(founders.R2))) founders.R2[is.nan(founders.R2)] <- 0
    
    # Estimate the founders MAF if given
    snp_maf <- unique(rbind(as.matrix(maf.df[,c("SNP_A", "MAF_A")]), as.matrix(maf.df[,c("SNP_B", "MAF_B")])))
    founders.maf <- snp_maf[which(snp_maf[,1] %in% founders[1]), 2]
    if(length(founders.maf) == 0) founders.maf <- 0.5
    
    
    # Generate non-homolgous markov models for the founder SNPs
    mcCCRC <- generateNHMC(founders, stateNames, founders.R2)
    
    # Calculates the probabilities of founder-SNP to sub-founder SNPs
    subf.p <- lapply(seq_along(founders.idx[-length(founders.idx)]), generateBinProb, founders.idx, ld.mat)
    
    # Extracts all SNP positions in the haplotype block between first and last founder SNP
    pos.df <- filterPos(pos.df, all.ids=colnames(ld.mat), fidx=founders.idx)
    
    mc.models <- list("subfounder"=subf.p,
                      "founder"=mcCCRC,
                      "fmaf"=founders.maf,
                      "ids"=pos.df$V3)  
    save(mc.models, file=mcmodel.id)
    write.table(pos.df, file=paste0(opt$id, "_filt.pos"), quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE)
    
    if(RUNHELPER){
      file.remove(paste0(opt$dir, 
                         paste0("/", opt$id, 
                                c(".log", ".maf.log", ".nosex", 
                                  ".tags.list", "_founders.txt", 
                                  ".ld.gz", ".maf.ld.gz", ".pos", ".maf.nosex"))))
    }
  }
  
  ## Simulations:
  if(grepl("sample", opt$mode, ignore.case=TRUE)){
    load(mcmodel.id)
    
    founders.maf <- mc.models[['fmaf']]
    subf.p <- mc.models[['subfounder']]
    mcCCRC <- mc.models[['founder']]
    ids <- mc.models[['ids']]
    
    subf.rand <- lapply(subf.p, sampleFromBinomial, opt$weight)
    
    m.sim <- sampleFromNHMC(mcCCRC, founders.maf, 100, stateNames)
    
    sims <- combineSampling(m.sim, stateNames, subf.rand, ids)
    
    write.table(sims, file=paste0(opt$id, "_simulations.txt"), sep="\t", col.names=TRUE, quote=FALSE)
  }
  
  if(!is.null(err)){
    cat(paste0("ERROR: Issues modelling haploblock ", each.row, "/", err, "...\n"))
  } else {
    cat(paste0("Finished modelling haploblock ", each.row, "/", max(iter.range), "...\n"))
  }
  err
}


#########
### MAIN

if(is.null(opt$dir)) opt$dir <- getwd()
stateNames = c(opt$ref, opt$alt)

setwd(opt$dir)

RUNHELPER <-FALSE
if(!is.null(opt$geno) && !is.null(opt$genoid)) RUNHELPER <- TRUE

if(RUNHELPER){
  base.dir <- gsub("/haploblock/.*", "", opt$dir)
  blocks.det <- file.path(base.dir, paste0(opt$genoid, ".blocks.det"))
  blocks.det <- read.table(blocks.det, header=TRUE, stringsAsFactors=FALSE, check.names=FALSE)
  blocks.det$range <- paste0(blocks.det$CHR, ":", blocks.det$BP1, "-", blocks.det$BP2)
  iter.range <- c(1:nrow(blocks.det))
} else {
  iter.range <- 1
}

err.row <- sapply(iter.range, function(each.row) LDMC(each.row))
err.row <- unlist(err.row)
if(any(!is.null(err.row))) sapply(err.row, function(each.row) LDMC(each.row))
  


