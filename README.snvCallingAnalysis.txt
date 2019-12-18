== Description

The SNV workflow takes two merged bam files as input (control and tumor sample)
It performs four steps:
1) snv calling
2) snv annotation
3) snv deep annotation
4) snv filter

== Run flags / switches

Switch                      Default Description
runSNVMetaCallingStep       false   Run a single job for snv calling instead of a job per chromosome.
                                    The meta job is optimized in regards to runtime and cpu / memory utilization.
                                    On the other hand it is a larger job which might be more difficult to schedule.
runDeepAnnotation           true    Run the deep annotation step or stop the workflow before it.
runFilter                   true    Run the filter step or stop the workflow before it.
runOnPancan                 false   Run a special analysis type for pancancer type projects.

== Changelist
* Version update to 1.2.166-3
  Added feature to the script snv_extractor.pl and the bash file filter_vcf.sh. Introduced a "Whitelist" output based on the 
  tab delimited file /icgc/ngs_share/assemblies/hg19_GRCh37_1000genomes/CustomGeneLists/INFORM_druggable_genes_ensembl_geneIDs_20191205.tsv
  Exonic somatic variants are filtered based on a more lenient confidence threshold that can be specified at the commandline (default value is 5)

* Version update to 1.2.166-2 

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
