# DKFZ SNVCalling Workflow

An SNV calling workflow developed in the Applied Bioinformatics and Theoretical Bioinformatics groups at the DKFZ. An earlier version (pre Github) of this workflow was used in the [Pancancer](https://dockstore.org/containers/quay.io/pancancer/pcawg-dkfz-workflow) project.

> <table><tr><td><a href="https://www.denbi.de/"><img src="docs/images/denbi.png" alt="de.NBI logo" width="300" align="left"></a></td><td><strong>Your opinion matters!</strong> The development of this workflow is supported by the <a href="https://www.denbi.de/">German Network for Bioinformatic Infrastructure (de.NBI)</a>. By completing <a href="https://www.surveymonkey.de/r/denbi-service?sc=hd-hub&tool=SNVCallingWorkflow">this very short (30-60 seconds) survey</a> you support our efforts to improve this tool.</td></tr></table>

## Installation

To run the workflow you first need to install a number of components and dependencies.

* You need a working [Roddy](https://github.com/TheRoddyWMS/Roddy) installation. The version depends on the workflow version you want to use. You can find it in the [buildinfo.txt](buildinfo.txt) under 'RoddyAPIVersion'. Please follow the instructions for the installation of Roddy itself and the PluginBase and the DefaultPlugin components. The main reference here is the [Roddy documentation](https://roddy-documentation.readthedocs.io/installationGuide.html).
* Install the version you need -- either from the release tarballs or with git clone into your plugin directory.

Furthermore you need a number of tools and of course reference data, like a genome assembly and annotation databases.

### Tool installation

The workflow contains a description of a [Conda](https://conda.io/docs/) environment. A number of Conda packages from [BioConda](https://bioconda.github.io/index.html) are required. 

First install the BioConda channels:
```
conda config --add channels r
conda config --add channels defaults
conda config --add channels conda-forge
conda config --add channels bioconda
conda config --add channels bioconda-legacy
```

Then install the environment

```
conda env create -n SNVCallingWorkflow -f $PATH_TO_PLUGIN_DIRECTORY/resources/analysisTools/snvPipeline/environments/conda.yml
```

The name of the Conda environment is arbitrary but needs to be consistent with the `condaEnvironmentName` variable. The default for that variable is set in `resources/configurationFiles/analysisSNVCalling.xml`.

Note that the Conda environment not exactly the same as the software stack used for the [Pancancer](https://dockstore.org/containers/quay.io/pancancer/pcawg-dkfz-workflow) project.

#### PyPy

PyPy is an alternative Python interpreter. Some of the Python scripts in the workflow can use PyPy to achieve higher performance by employing a fork of [hts-python](https://github.com/eilslabs/hts-python). Currently, this is not implemented for the Conda environment. For most cases you therefore should set the `PYPY_OR_PYTHON_BINARY` variable to just `python` to use the Python binary from the Conda environment. You could set up a `resources/analysisTools/snvPipeline/environments/conda_snvAnnotation.sh` similar to the `tbi-lsf-cluster_snvAnnotation.sh` file in the same directory.  

### Reference data installation

TBD

# Running the workflow

The SNV workflow takes two merged bam files as input (control and tumor sample)

It performs four steps:
   1. snv calling
   2. snv annotation
   3. snv deep annotation
   4. snv filter

## Configuration Values

|Switch                    |  Default     | Description
|--------------------------|--------------|-----------------------------------------------|
| bamfile_list             | empty        | Semicolon-separated list of BAM files, starting with the control's BAM. Each BAM file needs an index file with the same name as the BAM, but ".bai" suffixed |
| sample_list              | empty        | Semicolon-separated list of sample names in the same order as `bamfile_list` |
| possibleTumorSampleNamePrefixes | "( tumor )" | Bash-array of tumor sample name prefixes |
| possibleControlSampleNamePrefixes | "( control )" | Bash-array of control sample name prefixes |
| CHROMOSOME_INDICES | empty | Bash-array of chromosome names to which the analysis should be restricted |
| CHROMOSOME_LENGTH_FILE  | empty | Headerless TSV file with chromosome name, chromosome size columns |
| CHR_SUFFIX | "" | Suffix added to the chromosome names |
| CHR_PREFIX | "" | Prefix added to the chromosome names |
| extractSamplesFromOutputFiles | true | Refer to the documentation of the [COWorkflowBasePlugin](https://github.com/DKFZ-ODCF/COWorkflowsBasePlugin) for further information |
| PYPY_OR_PYTHON_BINARY | pypy | The binary to use for a some of the Python scripts. For `filter_PEoverlap.py` using a PyPy binary here also triggers the use of [hts-python](https://github.com/pjb7687/hts-python) instead of pysam.|  
| runSNVMetaCallingStep | false | Run a single job for snv calling instead of a job per chromosome. The meta job is optimized in regards to runtime and cpu / memory utilization. On the other hand it is a larger job which might be more difficult to schedule. |
| runDeepAnnotation | true | Run the deep annotation step or stop the workflow before it. |
| runFilter | true |Run the filter step or stop the workflow before it. |
| runOnPancan | false | Run a special analysis type for pancancer type  projects. | 
| NUMBER_OF_MISMATCHES_THRESHOLD | -1 | resources/analysisTools/snvPipeline/snvAnnotation.sh: Number of mismatches that are allowed per read in order to consider this read. |


## Example Call

```bash
roddy.sh run projectConfigurationName@analysisName patientId \
--useconfig=/path/to/your/applicationProperties.ini --configurationDirectories=/path/to/your/projectConfigs \
--useiodir=/input/directory,/output/directory/snv \
--usePluginVersion=SNVCallingWorkflow:1.3.2 \
--cvalues="bamfile_list:/path/to/your/control.bam;/path/to/your/tumor.bam,sample_list:normal;tumor,possibleTumorSampleNamePrefixes:tumor,possibleControlSampleNamePrefixes:normal,REFERENCE_GENOME:/reference/data/hs37d5_PhiX.fa,CHROMOSOME_LENGTH_FILE:/reference/data/hs37d5_PhiX.chromSizes,extractSamplesFromOutputFiles:false"
```

### No Control

TBD


### Cross-Species Contaminations

In coding regions, the expected proportion of synonymous mutations compared to the total number of mutations should be
low. By contrast, a high proportion of synonymous mutations suggests cross-species contamination. Any value
 above 0.5 (i.e. at least 50% of mutations are synonymous) is indicating a contamination. A value below 0.35 is considered to be
  OK. Values in the range of 0.35-0.5 are unclear. 


## Changelog

* 2.1.1

  * Patch level: Removed hardcoded path to ExAC database

* 2.1.0

* 2.0.0

  - Upgrading to BCFtools/SAMtools/htslib to 1.9
  - Calling variants using BCFtools mpileup and samtools mpileup for lookup in control (similar to the old way)
  - PE overlaps are removed as before
  - Removing supplementary reads and left-out duplicate reads from variant calling
  - Deactivating the use to ExAC for no-control workflow filtering.
  - Updating to new local control and change the max AF threshold to 0.01

* 2.0.0

* 1.4.1

* 1.4.0

* 1.3.2

* 1.3.1

* 1.3.0

* 1.2.166-3

  * Long-term support version.

* 1.2.166-3

* 1.2.166-2

* 1.2.166-1

* 1.1.4-2

* 1.1.4-1

* 1.1.4

* 1.1.3

* 1.1.2

* 1.0.166-2_R2.4

* 1.0.166-1_priorLsfTransition

* 1.0.166-1

The versions listed below are not reflected in this repository and are only listed as legacy versions. At the time the code was managed in an SVN repository and its history was not completely migrated into the current git repository.

* 1.0.166

* 1.0.164

* 1.0.163

* 1.0.162

* 1.0.158

* 1.0.156

* 1.0.155

  * Move SNV Calling workflow to a custom plugin

* 1.0.132

* 1.0.131

  * Change workflow class to override another execute method. This makes the workflow a bit cleaner.
  * filterVcfForBias.py
  * filter_PEoverlap.py

* 1.0.114

* 1.0.109

  - Bugfix: Rainfall plots work now also if there are less than two SNVs on chromosome 1 or any other chromosome

* 1.0.105

  - Bugfix: Use CHROMOSOME_LENGTH_FILE instead of CHROM_SIZES_FILE in snv filter (filter_vcf.sh) script.

* 1.0.104

* 1.0.103


## Contributors

Have a look at the [Contributors file](CONTRIBUTORS.md).
