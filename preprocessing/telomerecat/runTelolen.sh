#!/bin/bash
#
#$ -cwd
# module load telomerecat/3.1.2

PDIR=$1
mkdir -p ${PDIR}/output/telolength
cd ${PDIR}/output/telobam

echo "Converting telbam to length..."
telomerecat telbam2length \
-p 4 \
-v 2 \
-N 100 \
--output ${PDIR}/output/telolength/allTelbamLengths.csv \
$(ls . | grep "telbam.bam$")
