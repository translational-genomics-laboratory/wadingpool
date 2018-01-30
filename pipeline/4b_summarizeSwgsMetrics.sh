mkdir -p ${PDIR}/qc/picard/scripts ${PDIR}/qc/picard/output_summary
cd ${PDIR}/qc/picard/scripts

ln -s ${GIT}/analysis/qc/collapseAll.indel.sh .
ln -s ${GIT}/analysis/qc/collapseAll.wgs.sh .
ln -s ${GIT}/analysis/qc/plotPicardWgsMetrics.R .
ln -s ${GIT}/analysis/qc/runPlotPicard.sh .

echo "${bold}WadingPool:${normal} Sumarizing the QC metrics...";
sh collapseAll.indel.sh ${PDIR}"/qc/picard/"
sh collapseAll.wgs.sh ${PDIR}"/qc/picard/"
ln -s ${PDIR}"/qc/picard/output_iSize/allOL.isize.txt" ${PDIR}"/qc/picard/output_summary"
ln -s ${PDIR}"/qc/picard/output_wgs/allOL.wgsMetrics.txt" ${PDIR}"/qc/picard/output_summary"
sh runPlotPicard.sh ${PDIR}"/qc/picard/"
