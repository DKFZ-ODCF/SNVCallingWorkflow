# DKFZ SNVCalling Workflow

An SNV calling workflow developed in the Applied Bioinformatics and Theoretical Bioinformatics groups at the DKFZ. On earlier version (pre Github) of this workflow was used in the [Pancancer](https://github.com/TheRoddyWMS/BatchEuphoria/pull/124) project.

> <table><tr><td><a href="https://www.denbi.de/"><img src="docs/images/denbi.png" alt="de.NBI logo" width="300" align="left"></a></td><td><strong>Your opinion matters!</strong> The development of this workflow is supported by the <a href="https://www.denbi.de/">German Network for Bioinformatic Infrastructure (de.NBI)</a>. By completing <a href="yet unknown">this very short survey</a> you support our efforts to improve this tool.</td></tr></table>

## Installation

To run the workflow you first need to install a number of components and dependencies.

* You need a working [Roddy](https://github.com/TheRoddyWMS/Roddy) installation. The version depends on the workflow version you want to use. You can find it in the [buildinfo.txt](buildinfo.txt) under 'RoddyAPIVersion'. Please follow the instructions for the installation of Roddy itself and the PluginBase and the DefaultPlugin components. There reference here is the [Roddy documentation](https://roddy-documentation.readthedocs.io/installationGuide.html).
* Install the version you need -- either from the release tarballs or with git clone into your plugin directory.

Furthermore you need a number of tools and of course reference data, like a genome assembly and annotation databases.

### Tool installation

The workflow uses Conda as tool to provide the required tools. Note that the Conda environment YAML files tend to outdate, such that at some point some package is lost from all standard channels. We are working on providing a DKFZ Conda channel for long term provisioning of DKFZ workflows.  

#### PyPy

PyPy is an alternative Python interpreter. Some of the Python scripts in the workflow can use PyPy to achieve higher performance. Simply set the `PYPY_OR_PYTHON_BINARY` to a PyPy interpreter.

##### PyPy with DKFZ/LSF/Modules environment

Note that the following feature is currently not implemented in the Conda-based environment but only in the environment used at the DKFZ (based on LSF and environment modules). Further performance is gained through the employment of a PyPy/hts-python combination instead of a CPython/pysam in the `filter_PEoverlap.py` script. Again, simply set `PYPY_OR_PYTHON_BINARY` to a PyPy interpreter and the workflow uses the implementation based on hts-python.

### Reference data installation

Sebastians script

## Running the workflow

bamfile_list
Each BAM file needs to be accompanied by an index file with the same name but suffixed by .bai

### No Control

