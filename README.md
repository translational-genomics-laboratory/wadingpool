# WadingPool

WadingPool is a tool designed by the Pugh Lab (University Healthy Network, Toronto, ON) that is actively maintained and developed by Rene Quevedo.  This tool aims to create a pipeline facilitates the analysis of shallow whole-genome sequencing (< 0.3x coverage).  Included in this pipeline is tools to create copy-number profiles (via QDNAseq), purity estimation, LOH/zygosity inference, genotype identification and telomere status profiling.


## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.
Setting up the environment is best saved as a shell file labelled "0_setupEnvironment.sh", a sample file has been provided in the pipeline directory:
```
cat << EOF > 0_setupEnvironment.sh

GIT="/path/to/git/wadingpool"
PDIR='/path/to/working_directory'   # Directory you created for analysis
FASTQDIRS=${PDIR}'/fastq_dirs.txt'  # Pre-existing file containing absolute paths to directory harboring all raw fastq.gz files
IDLIST=${PDIR}'/id_list.txt'        # Generated in the "1_runSymlinkAlignFastq.sh" script
GENERATEIDS=true
LIBID='180124_NB501085_0191_AH373KBGX5_Project_1'    # Unique library ID
REFERENCEFILES='/mnt/work1/users/pughlab/references/igenome/coclean_reference_files.sh'
DBSNPDIR='/mnt/work1/users/pughlab/references/dbsnp/dbsnp_b146_GRCh37p13/common/bed'
DBSNPID='common_all_20151104.bed'

EOF
```

To facilitate the pre-processing steps, there are several helper bash-scripts that set-up the directory structure and execute the relevant commands.  Each helper script is aptly named to be descriptive of its function; more information of each script can be found as well.  Additionally, 
```
# Sets up all the environment variables requires for the pipeline
source 0_setupEnvironment.sh

# Processes the fastq > bam alignment;  Any pipeline can replace these parts, as long as the output is a chromosome-subsetted bam file
sh ${GIT}/pipeline/1_symlinkAndAlignFastq.sh
sh ${GIT}/pipeline/2_cocleanAllBams.sh
sh ${GIT}/pipeline/2b_symlinksCocleanedBams.sh
sh ${GIT}/pipeline/3_subsetByChromosomes.sh

# Runs all the QC steps to calculate coverage and general metrics
sh ${GIT}/pipeline/4_getSwgsMetrics.sh
sh ${GIT}/pipeline/4b_summarizeSwgsMetrics.sh

# Runs all the telomere-quantifying steps using TelomereCat
sh ${GIT}/pipeline/5_runTelomerecat.sh
sh ${GIT}/pipeline/5b_getTelomereLength.sh

# Runs all the variant calling steps, plus filtering for heterozygous or covered variants
sh ${GIT}/pipeline/6_variantCalling.sh
sh ${GIT}/pipeline/6b_selectHetVariants.sh
sh ${GIT}/pipeline/6c_mergeAllVcf.sh

# Runs all the variant calling steps, plus filtering for heterozygous or covered variants
sh ${GIT}/pipeline/6_variantCalling.sh
sh ${GIT}/pipeline/6b_selectHetVariants.sh
sh ${GIT}/pipeline/6c_mergeAllVcf.sh

# Runs QDNAseq and ichorCNA copy number calling on the bam files
sh ${GIT}/pipeline/7a_copyNumberCalling_QDNAseq.sh
sh ${GIT}/pipeline/7b_copyNumberCalling_ichorCNA.sh

```

### Prerequisites

All pre-processing steps are developed for the Sun Grid Engine

```
Give examples
```

### Installing

A step by step series of examples that tell you have to get a development env running

Say what the step will be

```
Give the example
```

And repeat

```
until finished
```

End with an example of getting some data out of the system or using it for a little demo

## Running the tests

Explain how to run the automated tests for this system

### Break down into end to end tests

Explain what these tests test and why

```
Give an example
```

### And coding style tests

Explain what these tests test and why

```
Give an example
```

## Deployment

Add additional notes about how to deploy this on a live system

## Built With

* [Dropwizard](http://www.dropwizard.io/1.0.2/docs/) - The web framework used
* [Maven](https://maven.apache.org/) - Dependency Management
* [ROME](https://rometools.github.io/rome/) - Used to generate RSS Feeds

## Contributing

Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.

## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://github.com/your/project/tags). 

## Authors

* **Rene Quevedo** - *Initial work* 

See also the list of [contributors](https://github.com/your/project/contributors) who participated in this project.

## License

This project is licensed under 

## Acknowledgments

* [QDNAseq](https://github.com/your/project/contributors)
* Inspiration
* etc

