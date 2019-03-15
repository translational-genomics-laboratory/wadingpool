module load R-gsi/3.5.1

pdir=$1

Rscript plotPicardWgsMetrics.R ${pdir}/output_summary/allOL.wgsMetrics.txt \
${pdir}/output_summary/allOL.isize.txt \
${pdir}/output_summary/swgs.cnMetrics.pdf
