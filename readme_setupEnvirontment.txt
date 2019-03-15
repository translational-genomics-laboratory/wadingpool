
cat << EOF > 0_setupEnvironment.sh
RUNID="180921_M00146_0089_000000000-D4WHY"
BINSIZE=50
REGEX="bam$"
BWAINDEX="/.mounts/labs/PDE/data/reference/hg19_random/fasta/UCSC/hg19_random.fa"
ICHORPATH="/.mounts/labs/TGL/gsi/tools/ichorCNA/scripts/runIchorCNA.R"
READCOPY="/.mounts/labs/TGL/gsi/tools/hmmcopy_utils/bin/readCounter"
GIT="/.mounts/labs/TGL/gsi/tools/wadingpool"
PDIR="/scratch2/groups/tgl/wadingPool"  # Directory you created for analysis
FASTQDIRS="${PDIR}/fastq_dirs.txt"  # Pre-existing file containing absolute paths to directory harboring all raw fastq.gz files
IDLIST="${PDIR}/id_list.txt"      # Generated in the "1_runSymlinkAlignFastq.sh" script
GENERATEIDS=true
LIBID="LIB26036"    # Unique library ID
REFERENCEFILES="/mnt/work1/users/pughlab/references/igenome/coclean_reference_files.sh"
DBSNPDIR="/mnt/work1/users/pughlab/references/dbsnp/dbsnp_b146_GRCh37p13/common/bed"
DBSNPID="common_all_20151104.bed"
EOF