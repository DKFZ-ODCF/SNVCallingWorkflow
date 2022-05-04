# Description

The SNV workflow takes two merged bam files as input (control and tumor sample)
It performs four steps:

1. SNV calling
2. SNV annotation
3. SNV deep annotation 
4. SNV filter

## Installation of the virtual env at central place in the cluster

This assumes, `tbiLsfVirtualEnvDir` is accessible from all compute nodes. Furthermore, Python 2.7.9 in this example is assumed to be available via a module named "python/2.7.9".

```bash
module load python/2.7.9
virtualenv "$tbiLsfVirtualEnvDir"
source "$tbiLsfVirtualEnvDir/bin/activate"
pip install -r "$pluginInstallationDir/resources/analysisTools/snvPipeline/environments/requirements.txt
```

Then, in your configuration files, you need to set the `tbiLsfVirtualEnvDir` variable.


## Run flags / switches

Switch                      Default Description
runSNVMetaCallingStep       false   Run a single job for snv calling instead of a job per chromosome.
                                    The meta job is optimized in regards to runtime and cpu / memory utilization.
                                    On the other hand it is a larger job which might be more difficult to schedule.
runDeepAnnotation           true    Run the deep annotation step or stop the workflow before it.
runFilter                   true    Run the filter step or stop the workflow before it.
runOnPancan                 false   Run a special analysis type for pancancer type projects.

## Changelist

* Version update to 1.2.166-4

  * introduced variable `tbiLsfVirtualEnvDir` to configure the path to the virtual env used in the `tbi-lsf-cluster.sh` environment setup script.

* Version update to 1.2.166-3

* Version update to 1.2.166-2

* Version update to 1.2.166-1

* Version update to 1.2.166

* Version update to 1.0.166

* Version update to 1.0.164

* Version update to 1.0.163

* Version update to 1.0.162

* Version update to 1.0.158

* Version update to 1.0.156

* Version update to 1.0.155

Move SNV Calling workflow to a custom plugin

* Version update to 1.0.132

* Version update to 1.0.131

- Change workflow class to override another execute method. This makes the workflow a bit cleaner.

- filterVcfForBias.py:

- filter_PEoverlap.py:

* Version update to 1.0.114

* Version update to 1.0.109

- Bugfix: Rainfall plots work now also if there are less than two SNVs on chromosome 1 or any other chromosome

* Version update to 1.0.105

- Bugfix: Use CHROMOSOME_LENGTH_FILE instead of CHROM_SIZES_FILE in snv filter (filter_vcf.sh) script.

* Version update to 1.0.104

* Version update to 1.0.103
