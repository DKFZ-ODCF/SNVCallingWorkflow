# DKFZ SNVCalling Workflow

An SNV calling workflow developed in the Applied Bioinformatics and Theoretical Bioinformatics groups at the DKFZ. An earlier version (pre Github) of this workflow was used in the [Pancancer](https://dockstore.org/containers/quay.io/pancancer/pcawg-dkfz-workflow) project.

> <table><tr><td><a href="https://www.denbi.de/"><img src="docs/images/denbi.png" alt="de.NBI logo" width="300" align="left"></a></td><td><strong>Your opinion matters!</strong> The development of this workflow is supported by the <a href="https://www.denbi.de/">German Network for Bioinformatic Infrastructure (de.NBI)</a>. By completing <a href="https://www.surveymonkey.de/r/denbi-service?sc=hd-hub&tool=SNVCallingWorkflow">this very short (30-60 seconds) survey</a> you support our efforts to improve this tool.</td></tr></table>

## Citing

The SNV workflow was in the pan-cancer analysis of whole genomes (PCAWG) and can be cited with the following publication:

* Pan-cancer analysis of whole genomes.<br>
  The ICGC/TCGA Pan-Cancer Analysis of Whole Genomes Consortium.<br>
  Nature volume 578, pages 82–93 (2020).<br>
  DOI [10.1038/s41586-020-1969-6](https://doi.org/10.1038/s41586-020-1969-6)

Containers are available in [Dockstore](https://dockstore.org/containers/quay.io/pancancer/pcawg-dkfz-workflow).

Furthermore, the workflow is regularly used at the DKFZ throug the automation system [One-Touch Pipeline](https://gitlab.com/one-touch-pipeline/otp):

* OTP: An automatized system for managing and processing NGS data.<br>
  Eva Reisinger, Lena Genthner, Jules Kerssemakers, Philip Kenschea, Stefan Borufka, Alke Jugold, Andreas Kling, Manuel Prinz, Ingrid Scholz, Gideon Zipprich, Roland Eils, Christian Lawerenz, Jürgen Eils.<br>
  Journal of Biotechnology, volume 261, pages 53-62 (2017).<br>
  DOI: [10.1016/j.jbiotec.2017.08.006](https://doi.org/10.1016/j.jbiotec.2017.08.006)

## Installation

To run the workflow you first need to install a number of components and dependencies.

* You need a working [Roddy](https://github.com/TheRoddyWMS/Roddy) installation. The version depends on the workflow version you want to use. You can find it in the [buildinfo.txt](buildinfo.txt) under 'RoddyAPIVersion'. Please follow the instructions for the installation of Roddy itself and the PluginBase and the DefaultPlugin components. The main reference here is the [Roddy documentation](https://roddy-documentation.readthedocs.io/installationGuide.html).
* Install the version you need -- either from the release tarballs or with git clone into your plugin directory.

Furthermore you need a number of tools and of course reference data, like a genome assembly and annotation databases.

### Tool installation

#### Conda

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

###### VirtualEnv

The virtual env does not provide a full software stack. This section is only relevant for the ODCF cluster and mostly for developers. In rare cases it may be necessary for other users to set the `tbiLsfVirtualEnvDir`.

This assumes, `tbiLsfVirtualEnvDir` is installed centrally and accessible from all compute nodes. Furthermore, Python 2.7.9 in this example is assumed to be available via a module named "python/2.7.9".

```bash
module load python/2.7.9
virtualenv "$tbiLsfVirtualEnvDir"
source "$tbiLsfVirtualEnvDir/bin/activate"
pip install -r "$pluginInstallationDir/resources/analysisTools/snvPipeline/environments/requirements.txt
```

Then, in your configuration files, you need to set the `tbiLsfVirtualEnvDir` variable.


#### PyPy

> PyPy support is broken and deprecated (2.2.0).

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

| Switch                            | Default       | Description                                                                                                                                                                                                                     |
|-----------------------------------|---------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| bamfile_list                      | empty         | Semicolon-separated list of BAM files, starting with the control's BAM. Each BAM file needs an index file with the same name as the BAM, but ".bai" suffixed                                                                    |
| sample_list                       | empty         | Semicolon-separated list of sample names in the same order as `bamfile_list`                                                                                                                                                    |
| possibleTumorSampleNamePrefixes   | "( tumor )"   | Bash-array of tumor sample name prefixes                                                                                                                                                                                        |
| possibleControlSampleNamePrefixes | "( control )" | Bash-array of control sample name prefixes                                                                                                                                                                                      |
| CHROMOSOME_INDICES                | empty         | Bash-array of chromosome names to which the analysis should be restricted                                                                                                                                                       |
| CHROMOSOME_LENGTH_FILE            | empty         | Headerless TSV file with chromosome name, chromosome size columns                                                                                                                                                               |
| CHR_SUFFIX                        | ""            | Suffix added to the chromosome names                                                                                                                                                                                            |
| CHR_PREFIX                        | ""            | Prefix added to the chromosome names                                                                                                                                                                                            |
| extractSamplesFromOutputFiles     | true          | Refer to the documentation of the [COWorkflowBasePlugin](https://github.com/DKFZ-ODCF/COWorkflowsBasePlugin) for further information                                                                                            |
| PYPY_OR_PYTHON_BINARY             | python        | Leave this set to "python". Some of the Python scripts used to work with PyPy (e.g. `filter_PEoverlap.py` that also uses [hts-python](https://github.com/pjb7687/hts-python) instead of pysam), but this is deprecated.         |
| runSNVMetaCallingStep             | false         | Run a single job for snv calling instead of a job per chromosome. The meta job is optimized in regards to runtime and cpu / memory utilization. On the other hand it is a larger job which might be more difficult to schedule. |
| runDeepAnnotation                 | true          | Run the deep annotation step or stop the workflow before it.                                                                                                                                                                    |
| runFilter                         | true          | Run the filter step or stop the workflow before it.                                                                                                                                                                             |
| runOnPancan                       | false         | Run a special analysis type for pancancer type projects. This will produce an additional VCF.                                                                                                                                   |
| NUMBER_OF_MISMATCHES_THRESHOLD    | -1            | resources/analysisTools/snvPipeline/snvAnnotation.sh: Number of mismatches that are allowed per read in order to consider this read.                                                                                            |

Please have a look at the `resources/configurationFiles/analysisSNVCalling.xml` for a more complete list of parameters.

## Example Call

```bash
roddy.sh run projectConfigurationName@analysisName patientId \
  --useconfig=/path/to/your/applicationProperties.ini \
  --configurationDirectories=/path/to/your/projectConfigs \
  --useiodir=/input/directory,/output/directory/snv \
  --usePluginVersion=SNVCallingWorkflow:1.3.2 \
  --cvalues="bamfile_list:/path/to/your/control.bam;/path/to/your/tumor.bam,sample_list:normal;tumor,possibleTumorSampleNamePrefixes:tumor,possibleControlSampleNamePrefixes:normal,REFERENCE_GENOME:/reference/data/hs37d5_PhiX.fa,CHROMOSOME_LENGTH_FILE:/reference/data/hs37d5_PhiX.chromSizes,extractSamplesFromOutputFiles:false"
```

### No Control

* Set the parameter `isNoControlWorkflow` to `true`.
* In `bamfile_list` provide only the tumor BAM as only BAM file. The same for `sample_list`.
* `possibleTumorSampleNamePrefixes` must still be set to the name of the tumor sample

### Cross-Species Contamination

In coding regions, the expected proportion of synonymous mutations compared to the total number of mutations should be low. By contrast, a high proportion of synonymous mutations suggests cross-species contamination. Any value above 0.5 (i.e. at least 50% of mutations are synonymous) is indicating a contamination. A value below 0.35 is considered to be OK. Values in the range of 0.35-0.5 are unclear. 


### VCF Conversion Script (branch: ReleaseBranch_1.2.166)

The [convertToStdVCF.py](./blob/master/resources/analysisTools/snvPipeline/convertToStdVCF.py) may be helpful to convert the VCFs produced by this workflow into standard 4.2-conform VCFs.

There are, however, few caveats you should be aware of:

  * This is basically rescued old code that never has been extensively tested for 100% conformance in all cases. We advise you check the resulting VCFs with a [VCF validator](https://github.com/EBIvariation/vcf-validator).
  * Not all columns of the DKFZ VCF are currently migrated into the standard-conform VCF. Some columns may be lost.
  * The VCF must only contain a single sample column.
  * Only uncompressed VCFs are processed and the input must not be a stream/pipe, because the script does two passes over the input.

```bash
convertToStdVCF.py <infile> <outfile> <sampleId> [<config.json>]
```

The optional configuration JSON file defaults to the `convertToStdVCF.json` residing besides the `convertToStdVCF.py` script. The [default JSON](./blob/master/resources/analysisTools/snvPipeline/convertToStdVCF.json) contains some comments that may help you to customize the script output.

## Changelog

* upcoming

  * patch: Remove all code related to PyPy and hts-python (including `copysam.py` and `PYPY_OR_PYTHON_BINARY`)

* 2.2.0

  * minor: Update virtualenv
    * Updated environment to matplotlib 1.5.3. Version 1.4.3 seems to be incompatible to numpy 1.11.3 (now; `import matplotlib.pyplot` failed).
    * Used pysam 0.16.0.1 from the orginal 2.1.1 environment. This breaks the PyPy support (`copysam.py` is not yet adapted) but is necessary for the changes done for workflow version 2.1.1.
    * Updated to BioPython 1.71. This was lost in workflow version 2.1.1 but is necessary for some data `from bio import bgzf` is needed (`vcfparser.py`).
  * minor: PyPy is deprecated. Configuration values related to it are (mostly) removed.
  * minor: Allow configuration of virtualenv in `tbi-lsf-cluster.sh`
  * Lifted readme from ReleaseBranch_1.2.166
  * > Note: Again the Conda environment is broken and cannot get rebuild. 

* 2.1.1

  * patch: Removed hardcoded path to ExAC database

* 2.1.0

* 2.0.0

  - Upgrading to BCFtools/SAMtools/htslib to 1.9
  - Calling variants using BCFtools mpileup and samtools mpileup for lookup in control (similar to the old way)
  - PE overlaps are removed as before
  - Removing supplementary reads and left-out duplicate reads from variant calling
  - Deactivating the use to ExAC for no-control workflow filtering.
  - Updating to new local control and change the max AF threshold to 0.01

* 1.4.1

* 1.4.0

* 1.3.2

* 1.3.1

* 1.3.0

* 1.2.166-3 (ReleaseBranch_1.2.166)

  * Long-term support version.

* 1.2.166-3

* 1.2.166-2

* 1.2.166-1

  * This version derives from 1.0.166-2_R2.4 and was done in the context of a migration of our cluster from PBS to IBM LSF and associated update of the workflow manager Roddy.

* 1.1.4-2

* 1.1.4-1

* 1.1.4

* 1.1.3

* 1.1.2

* 1.0.166-2_R2.4

* 1.0.166-1_priorLsfTransition

* 1.0.166-1

The versions listed below are absent from this repository and are only listed as legacy versions. At the time, the code was managed in an SVN repository and its history was not migrated into the current git repository.

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
