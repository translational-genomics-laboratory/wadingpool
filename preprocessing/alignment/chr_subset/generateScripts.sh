#!/bin/bash
#
#$ -cwd
#$ -S /bin/bash
#
############################################################
bold=$(tput bold)
normal=$(tput sgr0)

## Specify PDIR and LIBID as passed in or hardcoded
if [ $# -eq 0 ]; then
  echo "${bold}WadingPool:${normal} Using hardcoded PDIR"
  PDIR='/mnt/work1/users/pughlab/projects/test_dir/swgs'
else
  echo "${bold}WadingPool:${normal} Using PDIR specified by configuration file"
  PDIR=$1
fi

rm queueJobs.sh
for id in $(cat id_list.txt); do
  for chrom in $(seq 1 22); do
  
    cat << EOF > subset.chr${chrom}.${id}.sh
#!/bin/bash
#
#$ -cwd

# . /etc/bashrc
module load samtools/1.2 

chr="chr${chrom}"
PDIR=${PDIR}
BAMID="${id}.bam"

inDir="\${PDIR}/bam"
inBam="\${inDir}/\${BAMID}"

outDir="\${PDIR}/chr_subset/\${chr}"
mkdir -p \${outDir} 
cd \$outDir

echo \${chr}
samtools view -Sbh  \${inBam} "\${chr}" > \${chr}.\${BAMID}
samtools index \${chr}.\${BAMID}

EOF
  
    echo "qsub subset.chr${chrom}.${id}.sh" >> queueJobs.sh
  done
done