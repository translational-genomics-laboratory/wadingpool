## symlink FASTQ directory and align
bold=$(tput bold)
normal=$(tput sgr0)

cd ${PDIR}/data/cocleaned_bam
PDIR=${PDIR}/data/cocleaned_bam

mkdir symlinks
for i in $(ls .); do
  ln -s ${PDIR}/${i}/*ba[mi] ${PDIR}/symlinks;
done
