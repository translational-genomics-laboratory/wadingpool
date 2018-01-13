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
  PDIR='/mnt/work1/users/pughlab/projects/CHX/swgs'
else
  echo "${bold}WadingPool:${normal} Using PDIR specified by configuration file"
  PDIR=$1
fi
PDIR=${PDIR}'/telomerecat'

rm queueJobs.sh
for eachBam in $(cat id_list.txt); do
  cat << EOF > telocat.${eachBam}.sh
#!/bin/bash
#
#
#$ -cwd
#$ -S /bin/bash
module load telomerecat/3.1.2

PDIR=${PDIR}
mkdir \${PDIR}/output/telobam

echo "Making telbam..."
telomerecat bam2telbam \\
-v 2 \\
-p 4 \\
\${PDIR}/input/${eachBam}.bam

mv \${PDIR}/sh_scripts/${eachBam}_telbam.bam \${PDIR}/output/telobam

EOF
  echo "qsub telocat.${eachBam}.sh" >> queueJobs.sh
done