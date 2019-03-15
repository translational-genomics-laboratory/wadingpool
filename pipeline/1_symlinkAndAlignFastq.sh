## symlink FASTQ directory and align
bold=$(tput bold)
normal=$(tput sgr0)
cat << EOF
${bold}WadingPool:${normal} Symlink and Align FASTQs
    Helper script to take in the FASTQ directories specified by \$FASTQDIRS, symlink all the fastq.gz files, concatenate multiplexed FASTQ files, and then run bwa aligner on the files
    Assumes the fastq.gz files take the following identification structure:

      CHX-022-T_S44_L001_R1_001.fastq.gz
      [ID]_[LANE]_[READ]_[MULTIPLEX_ID].fastq.gz
      ID:CHX-022-T_S44
EOF


mkdir -p ${PDIR}/data
cd ${PDIR}/data
mkdir -p bam  chr_subset  cocleaned_bam  fastq  fastq_symlinks  logs  scripts  sh_scripts/bwa sh_scripts/coclean  symlink_dir

for i in $(cat $FASTQDIRS); do
  ln -s $i symlink_dir
  ln -s $i/*fastq.gz fastq_symlinks
  rm fastq_symlinks/\*fastq.gz
done

if($GENERATEIDS); then
  echo "${bold}WadingPool:${normal} A new id_list.txt will be generated based on the fastq's located in your fastq_dirs.txt directories";
  if [ -f ${IDLIST} ]; then
    mv ${IDLIST} ${IDLIST}.tmp
    echo "${bold}WadingPool:${normal} An existing id_list.txt has been found and been renamed as id_list.txt.tmp to avoid being overwritten, please rename appropriately.";
  fi
  echo "" > ${IDLIST}
  for i in $(cat $FASTQDIRS); do
    ls -1 $i | grep "fastq.gz" | sed "s/_L0.*//" | sort | uniq >> ${IDLIST}
  done
fi

cd sh_scripts/bwa
rm id_list.txt generateScripts.sh
ln -s ${IDLIST} .
ln -s ${GIT}/preprocessing/alignment/bwa/generateScripts.sh .
echo "${bold}WadingPool:${normal} Aligning fastq files...";
sh generateScripts.sh ${PDIR} ${LIBID} ${BWAINDEX}
sh queueJobs.sh
