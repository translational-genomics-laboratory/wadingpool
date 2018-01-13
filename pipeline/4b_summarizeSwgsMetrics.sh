mkdir -p ${PDIR}/qc/picard/scripts
cd ${PDIR}/qc/picard/scripts

ln -s ${GIT}/analysis/qc/collapseAll.indel.sh .
ln -s ${GIT}/analysis/qc/collapseAll.wgs.sh .
ln -s ${GIT}/analysis/qc/plotPicardWgsMetrics.R .
ln -s ${GIT}/analysis/qc/runPlotPicard.sh .

echo "${bold}WadingPool:${normal} Sumarizing the QC metrics...";
sh collapseAll.indel.sh ${PDIR}
sh collapseAll.wgs.sh ${PDIR}
sh runPlotPicard.sh ${PDIR}