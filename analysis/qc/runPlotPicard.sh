module load R/3.2.2

pdir=$1

Rscript plotPicardWgsMetrics.R ${pdir}/output_summary/allOL.wgsMetrics.txt \
${pdir}/output_summary/allOL.isize.txt \
${pdir}/output_summary/CHX.batch1.pdf