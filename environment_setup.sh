##############################
gatk_dir=/oicr/local/analysis/sw/gatk/GenomeAnalysisTK-3.5-0
JAVA=/oicr/local/lib/jvm/jdk1.6.0_25/bin/java
java=/.mounts/labs/PDE/Modules/sw/jvm/jdk1.8.0_91/bin/java
picard_jar=/.mounts/labs/PDE/Modules/sw/picard/2.12.1/picard.jar

##############################
REF=/.mounts/labs/PDE/data/gatkAnnotationResources/hg19_random.fa
grch371000gIndels=/.mounts/labs/TGL/gsi/databases/1000G_phase1.indels.b37.vcf_chr.vcf
grch37MillsIndels=/.mounts/labs/TGL/gsi/databases/Mills_and_1000G_gold_standard.indels.b37.vcf_chr.vcf
dbsnpVcf=/oicr/data/genomes/homo_sapiens_mc/dbSNP/hg19_random/Genomic/dbSNP147/dbsnp147_chr.vcf
BWAINDEX=/oicr/data/reference/genomes/homo_sapiens_mc/UCSC/hg19_random/Genomic/bwa/0.7.9/hg19_random.fa
###############################
SAMTOOLS=`which samtools`
###############################
#USER CUSTOM#
################################
RUNID="180921_M00146_0089_000000000-D4WHY"
BINSIZE=50
REGEX="bam$"
ICHORPATH="/.mounts/labs/TGL/gsi/tools/ichorCNA/scripts/runIchorCNA.R"
READCOPY="/.mounts/labs/TGL/gsi/tools/hmmcopy_utils/bin/readCounter"
GIT="/.mounts/labs/TGL/gsi/tools/WadingPool_Prisni"
PDIR="/scratch2/groups/tgl/wadingPool"  # Directory you created for analysis
FASTQDIRS="/scratch2/groups/tgl/wadingPool/fastq_dirs.txt"  # Pre-existing file containing absolute paths to directory harboring all raw fastq.gz files
PROV2="/.mounts/labs/seqprodbio/private/backups/seqware_files_report_latest.tsv.gz" #only for OICR
zgrep "TGL07" $PROV2 | grep "CASAVA" | grep "WG" | cut -f47 > $FASTQDIRS # only for OICR
IDLIST="/scratch2/groups/tgl/wadingPool/id_list.txt"      # Generated in the "1_runSymlinkAlignFastq.sh" script
GENERATEIDS=true
LIBID="LIB26036"    # Unique library ID
# DBSNPDIR="/mnt/work1/users/pughlab/references/dbsnp/dbsnp_b146_GRCh37p13/common/bed"
# DBSNPID="common_all_20151104.bed"
export PATH=/u/prath/anaconda2/bin:$PATH

R_LIBS=/.mounts/labs/gsiprojects/gsi/tools/R_libs


#####################################
# create synlinks of fastq files
FQSYM=${PDIR}/fastq_symlinks/
mkdir -p $FQSYM

BAMDIR=${PDIR}/bam/
mkdir -p $BAMDIR

LOGDIR=${PDIR}/logs
mkdir -p $LOGDIR

QC="$PDIR/QC"
mkdir -p $QC

PCQCDIR=$PDIR/Picard_metrics
mkdir -p $PCQCDIR

TLDIR=$PDIR/telomerecat
mkdir -p $TLDIR

TL=$TLDIR/telomere_length
mkdir -p $TL
