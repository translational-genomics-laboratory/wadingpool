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
PDIR=${PDIR}'/qc/picard'

rm queueJobs.sh
for id in $(cat id_list.txt); do
  cat << EOF > picardMetrics.$id.sh
#!/bin/bash
#
#$ -cwd

module load picard/2.4.1
module load igenome-human/hg19 

PDIR=${PDIR}
FILE="$id"

echo "Running CollectWgsMetrics..."
java -Xmx6g -jar \$picard_dir/picard.jar CollectWgsMetrics \\
I=\$PDIR/input/\$FILE.cocleaned.bam \\
O=\$PDIR/output_wgs/\$FILE.wgsMetrics.txt \\
R=\$REF


echo "Running insertSize Metrics..."
java -Xmx6g -jar \$picard_dir/picard.jar CollectInsertSizeMetrics \\
I=\$PDIR/input/\$FILE.cocleaned.bam \\
O=\$PDIR/output_iSize/\$FILE.isize.txt \\
H=\$PDIR/output_iSize/\$FILE.histogram.pdf \\
M=0.5

EOF

  echo "qsub picardMetrics."$id".sh" >> queueJobs.sh
done