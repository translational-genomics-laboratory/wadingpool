## symlink FASTQ directory and align
bold=$(tput bold)
normal=$(tput sgr0)

cd ${PDIR}/data/cocleaned_bam

mkdir symlinks
for i in $(ls .); do
  ln -s ${PDIR}/data/cocleaned_bam/${i}/*ba[mi] ${PDIR}/data/cocleaned_bam/symlinks;
done
rm ${PDIR}/data/cocleaned_bam/symlinks/\*ba\[mi\]