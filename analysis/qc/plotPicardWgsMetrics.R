args <- commandArgs(TRUE)
metricsFile <- args[1]
isizeFile <- args[2]
outPdf <- args[3]



isize.df <- read.table(isizeFile, sep="\t", header=TRUE, 
                       check.names=FALSE, stringsAsFactors=FALSE)
rownames(isize.df) <- gsub(".wgsMetrics.txt", "", isize.df$SAMPLES)
isize.df$SAMPLES <- NULL

hsmet.df <- read.table(metricsFile, sep="\t", header=TRUE, 
                       check.names=FALSE, stringsAsFactors=FALSE)
rownames(hsmet.df) <- gsub(".wgsMetrics.txt", "", hsmet.df$SAMPLES)
hsmet.df$SAMPLES <- NULL

# Plotting Features:
pdf(outPdf, height=12, width=10)
split.screen(c(3,1))
# Mean Coverage
screen(1)
par(mar=c(1,4.1, 2, 2.1))
cov.col.2 <- c(bait="gray")
if(max(hsmet.df[,c('MEAN_COVERAGE')]) < 0.01) yrange <- c(0, 0.01) else c(0, max(hsmet.df[,c('MEAN_COVERAGE')]))
barplot(t(as.matrix(hsmet.df[,c('MEAN_COVERAGE'), drop=FALSE])), 
        ylab="Mean coverage", xaxt='n', las=2, cex.names=0.75,
        ylim=yrange, col=cov.col.2)

screen(2)
par(mar=c(1,4.1, 2, 2.1))
barplot(t(as.matrix(isize.df[,c('MEDIAN_INSERT_SIZE'), drop=FALSE])), 
        ylab="Median insert size", xaxt='n', las=2, cex.names=0.75,
        ylim=c(0, 500), col=cov.col.2)

# Failed/filtered bases:
screen(3)
par(mar=c(6.1,4.1, 1, 2.1))
filt.col.5 <- c("#7fc97f", "#beaed4", "#fdc086", "#ffff99", "#386cb0")
barplot(t(hsmet.df[,c("PCT_EXC_MAPQ", "PCT_EXC_DUPE", "PCT_EXC_UNPAIRED",
                      "PCT_EXC_BASEQ", "PCT_EXC_OVERLAP")]),
        col=filt.col.5, cex.names=0.5, ylim=c(0, 1), las=2, ylab="Fraction of filtered aligned bases")
legend("topright", fill=filt.col.5, cex=1,
       legend=c("MAPQ < 20", "Duplicates", "Singletons", "BASEQ < 20", "Overlapping reads"))
close.screen(all.screens=TRUE)





# Percent of bases meeting coverage
targ.base.df <- hsmet.df[,c("PCT_1X", "PCT_5X", "PCT_10X",
                            "PCT_20X", "PCT_30X", "PCT_40X",
                            "PCT_50X", "PCT_100X")]
plot(0, type='n', ylim=c(0,0.5), xlim=c(0,8),
     ylab="Fraction of bases above coverage", xlab="Coverage", xaxt='n')
axis(1, at=c(1:dim(targ.base.df)[2]), labels=gsub("PCT_", "", colnames(targ.base.df)),
     las=2)
sample.col <- rainbow(dim(targ.base.df)[1])
for(each.row in c(1:dim(targ.base.df)[1])){
  lines(c(1:dim(targ.base.df)[2]),
        targ.base.df[each.row,], col=sample.col[each.row])
}
legend("topright", legend = rownames(targ.base.df), 
       lty=rep(1, dim(targ.base.df)[1]), lwd=rep(2.5, dim(targ.base.df)[1]),
       col = sample.col, cex=1)
dev.off()