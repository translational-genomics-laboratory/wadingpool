library(optparse)
option_list=list(
  make_option(c("-i", "--input"), type='character', default=NULL,
              help='Input filename; should be the output of Varscan/2.4.2 readcounts', metavar='character'),
  make_option(c("-r", "--reference"), type='character', default=NULL,
              help='Path to the matching detailed position file containing marker ID and ref/alt genotypes', metavar='character'),
  make_option(c("-o", "--output"), type='character', default='NULL',
              help='Output filename', metavar='character'),
  make_option(c("-p", "--path"), type='character', default='NULL',
              help='Path to the directory containing the input file and where to write the output', metavar='character')
)

opt_parser=OptionParser(usage = "Rscript %prog --input <INPUT FILE>
                        --output <OUTPUT FILE>
                        --path <OUTPUT DIRECTORY>", option_list=option_list)
opt=parse_args(opt_parser)

if(!is.null(opt$path)) setwd(opt$path)
det.pos <- opt$reference
rc.file <- opt$input
out.file <- opt$output

# Loading the data
cat("Loading in files...\n")
det.pos <- read.table(det.pos, header=TRUE, stringsAsFactors=FALSE, check.names=FALSE)
rc.df <- read.table(rc.file, header=TRUE, stringsAsFactors=FALSE,
                    colClasses = c(rep("character", 3), rep("integer", 2), rep("character", 2)),
                    check.names=FALSE, fill=TRUE)
colnames(rc.df) <- c("chrom", "position", "ref", "depth", "q20", "REF_INFO", "ALT_INFO")

# ALT info
cat("Getting Alt allele...\n")
rc.df$alt <- gsub(":.*", "", rc.df$ALT_INFO)
rc.df[which(rc.df$alt == ''),]$alt <- NA

# Populate depth info
cat("Getting depth info...\n")
rc.df$ref_depth <- as.integer(gsub("^.*?:(.*?):.*$", "\\1", rc.df$REF_INFO))
rc.df$alt_depth <- as.integer(gsub("^.*?:(.*?):.*$", "\\1", rc.df$ALT_INFO))
rc.df$alt_depth[is.na(rc.df$alt_depth)] <- 0

# Populate QUAL info
cat("Getting QUAL info...\n")
rc.df$ref_qual <- as.integer(gsub("^(.*?:){4}(.*?):.*$", "\\2", rc.df$REF_INFO, perl=TRUE)) 
rc.df$alt_qual <- as.integer(gsub("^(.*?:){4}(.*?):.*$", "\\2", rc.df$ALT_INFO, perl=TRUE)) 

# Add in marker ID
cat("Retrieving the SNP marker ID...\n")
det.pos <- as.data.frame(det.pos)
if(ncol(det.pos) != 5) stop("Positions dataframe assumes CHR POS ID REF ALT columns/ordering")
colnames(det.pos) <- c("chr", "pos", "id", "ref_det", "alt_det")
det.pos$uid <- gsub("^chr", "", paste0(det.pos[,1], "-", det.pos[,2]))
rc.df$uid <- gsub("^chr", "", paste0(rc.df$chrom, "-", rc.df$position))

# Merge the marker IDs and reform the dataframe
cat("Merging and writing...\n")
rc.df <- merge(rc.df, det.pos, by="uid", all.x = TRUE)
rc.adj.df <- rc.df[,c("chrom", "position", "id",
                      "ref", "alt", "ref_det", "alt_det",
                      "ref_depth", "alt_depth", "depth",
                      "ref_qual", "alt_qual")]    

# Write the updated table out
write.table(rc.adj.df, file=out.file, sep="\t",
            quote=FALSE, col.names=TRUE, row.names=FALSE)   

