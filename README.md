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

In coding regions, the expected proportion of synonymous mutations compared to the total number of mutations should be low. On the other hand, a high proportion of synonymous mutations may be an indicator of cross-species contamination. Any value above 0.5 (50% of mutations are synonymous) is a clear indicator of contamination. Anything below 0.35 is considered to be OK. Values in the range of 0.35-0.5 are unclear. 

## Contributors

Have a look at the [Contributors file](CONTRIBUTORS.md).
