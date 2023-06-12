# Convert DKFZ VCFs to standard-conform VCFs

The VCFs produced by the SNVCallingWorkflow are not standard conform in that some values are not added as additional columns after a single variant column. By contrast, in the [standard](https://samtools.github.io/hts-specs/) format, additional columns should only be used to show variants occurring in additional samples.

The `convertToStdVCF.py` script can be used to convert the DKFZ VCFs to standard VCFs (version 4.2).

## Execution

```bash
$ cat dkfz.vcf
##fileformat=VCFv4.1
##samtoolsVersion=0.1.19-44428cd
##reference=file:///icgc/ngs_share/assemblies/hg19_GRCh37_1000genomes/sequence/1KGRef/hs37d5.fa
##contig=<ID=1,length=249250621>
##INFO=<ID=DP,Number=1,Type=Integer,Description="Raw read depth">
##INFO=<ID=DP4,Number=4,Type=Integer,Description="# high-quality ref-forward bases, ref-reverse, alt-forward and alt-reverse bases">
##INFO=<ID=MQ,Number=1,Type=Integer,Description="Root-mean-square mapping quality of covering reads">
##INFO=<ID=FQ,Number=1,Type=Float,Description="Phred probability of all samples being the same">
##INFO=<ID=AF1,Number=1,Type=Float,Description="Max-likelihood estimate of the first ALT allele frequency (assuming HWE)">
##INFO=<ID=AC1,Number=1,Type=Float,Description="Max-likelihood estimate of the first ALT allele count (no HWE assumption)">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
##INFO=<ID=IS,Number=2,Type=Float,Description="Maximum number of reads supporting an indel and fraction of indel reads">
##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes for each ALT allele, in the same order as listed">
##INFO=<ID=G3,Number=3,Type=Float,Description="ML estimate of genotype frequencies">
##INFO=<ID=HWE,Number=1,Type=Float,Description="Chi^2 based HWE test P-value based on G3">
##INFO=<ID=CLR,Number=1,Type=Integer,Description="Log ratio of genotype likelihoods with and without the constraint">
##INFO=<ID=UGT,Number=1,Type=String,Description="The most probable unconstrained genotype configuration in the trio">
##INFO=<ID=CGT,Number=1,Type=String,Description="The most probable constrained genotype configuration in the trio">
##INFO=<ID=PV4,Number=4,Type=Float,Description="P-values for strand bias, baseQ bias, mapQ bias and tail distance bias">
##INFO=<ID=INDEL,Number=0,Type=Flag,Description="Indicates that the variant is an INDEL.">
##INFO=<ID=PC2,Number=2,Type=Integer,Description="Phred probability of the nonRef allele frequency in group1 samples being larger (,smaller) than in group2.">
##INFO=<ID=PCHI2,Number=1,Type=Float,Description="Posterior weighted chi^2 P-value for testing the association between group1 and group2 samples.">
##INFO=<ID=QCHI2,Number=1,Type=Integer,Description="Phred scaled PCHI2.">
##INFO=<ID=PR,Number=1,Type=Integer,Description="# permutations yielding a smaller PCHI2.">
##INFO=<ID=QBD,Number=1,Type=Float,Description="Quality by Depth: QUAL/#reads">
##INFO=<ID=RPB,Number=1,Type=Float,Description="Read Position Bias">
##INFO=<ID=MDV,Number=1,Type=Integer,Description="Maximum number of high-quality nonRef reads in samples">
##INFO=<ID=VDB,Number=1,Type=Float,Description="Variant Distance Bias (v2) for filtering splice-site artefacts in RNA-seq data. Note: this version may be broken.">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=GL,Number=3,Type=Float,Description="Likelihoods for RR,RA,AA genotypes (R=ref,A=alt)">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="# high-quality bases">
##FORMAT=<ID=DV,Number=1,Type=Integer,Description="# high-quality non-reference bases">
##FORMAT=<ID=SP,Number=1,Type=Integer,Description="Phred-scaled strand bias P-value">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="List of Phred-scaled genotype likelihoods">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	patient007	SEQUENCE_CONTEXT	INFO_control(VAF=variant_allele_fraction;TSR=total_variant_supporting_reads_incl_lowqual)	ANNOTATION_control	DBSNP	1K_GENOMES	ANNOVAR_FUNCTION	GENE	EXONIC_CLASSIFICATION	ANNOVAR_TRANSCRIPTS	SEGDUP	CYTOBAND	REPEAT_MASKER	DAC_BLACKLIST	DUKE_EXCLUDED	HISEQDEPTH	SELFCHAIN	MAPABILITY	SIMPLE_TANDEMREPEATS	CONFIDENCE	RECLASSIFICATION	seqBiasPresent_1	seqingBiasPresent_1	seqBiasPresent_2	seqingBiasPresent_2	Enhancers	CpGislands	TFBScons	ENCODE_DNASE	miRNAs_snoRNAs	miRBase18	COSMIC	miRNAtargets	CgiMountains	phastConsElem20bp	ENCODE_TFBS
1	10051	.	A	C	-0	RE;TAC;HSDEPTH;MAP;FRQ	DP=121;VDB=3.572478e-02;RPB=-7.344126e-01;AF1=0;AC1=0;DP4=72,5,1,2;MQ=41;FQ=-218;PV4=0.019,0.00027,1,0.11	GT:PL:GQ	1/1:0,191,255:99	CCTAACCCTA,CCCTAACCCT	DP=387;DP5=239,15,2,0,16;DP5all=269,75,2,3,38;ACGTNacgtnHQ=239,5,2,0,0,15,0,0,0,0;ACGTNacgtn=269,5,2,1,0,75,18,3,2,0;VAF=0.01;TSR=5;PBINOM=4.89276584292767e-78	unclear	.	.	intergenic	NONE(dist=NONE),DDX11L1(dist=1818)	.	.	Score=0.991956;Name=chr15:102446355	1p36.33	Simple_repeat_Simple_repeat_(CCCTAA)n	.	.	HSDR_1	4031;normScore=89.7	0.333333	trf_1;score=789	1	LQVSIG	.	.	.	.	.	.	.	.	.	.	.	.	CGI(10%);score=90	.	.
$ python2 convertToStdVCF.py -i dkfz.vcf -s test_sample
##fileformat=VCFv4.2
##INFO=<ID=ACGTNacgtnHQ_ctrl,Number=10,Type=Integer,Description="Nucleotide counts on ?">
##INFO=<ID=ACGTNacgtn_ctrl,Number=10,Type=Integer,Description="Nucleotide counts on ?">
##INFO=<ID=ANNOVARTR,Number=1,Type=String,Description="Details of non-synonymous variation's impact on protein">
##INFO=<ID=CONF,Number=1,Type=Integer,Description="Empirical confidence level of call, ten levels (1-10), >7=trustworthy">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total depth">
##INFO=<ID=DP4,Number=4,Type=Integer,Description="Counts of REF and ALT reads on forward and reverse strand">
##INFO=<ID=DP5_ctrl,Number=5,Type=Integer,Description="Counts of ?">
##INFO=<ID=DP_ctrl,Number=1,Type=Integer,Description="Total depth of control">
##INFO=<ID=EXONICCL,Number=1,Type=String,Description="Functional implication if variation is located with respect to exon">
##INFO=<ID=GENE,Number=1,Type=String,Description="GeneSymbol if variation overlaps gene">
##INFO=<ID=MQ,Number=1,Type=Float,Description="RMS mapping quality">
##INFO=<ID=RECL,Number=1,Type=String,Description="Reclassification of the germline/somatic state based on additional information">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype, /:unphased, |:phased">
##FORMAT=<ID=PL,Number=0,Type=Integer,Description="Phred-scaled genotype likelihoods rounded to closest integer">
##contig=<ID=1,length=249250621>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	patient007
1	10051	.	A	C	-0	RE;TAC;HSDEPTH;MAP;FRQ	DP4=72,5,1,2;ACGTNacgtnHQ_ctrl=239,5,2,0,0,15,0,0,0,0;RECL=LQVSIG;DP_ctrl=387;DP5_ctrl=239,15,2,0,16;ANNOVARTR=.;MQ=41;EXONICCL=.;CONF=1;GENE=NONE(dist=NONE),DDX11L1(dist=1818);ACGTNacgtn_ctrl=269,5,2,1,0,75,18,3,2,0;DP=121	GT:PL:GQ	1/1:0,191,255:99
```

## Configuration

The file `convertToStdVCF.json` specifies how the non-standard columns are converted into attributes.

In this file, entries starting with "__" are considered as comment and ignored by the script.

The file contains three sections:

  * "FILTERS": Input columns mappet to key/value fields in the "filters" column.
  * "FORMAT": Input columns mapped to key/value fields in the "format" column.
  * "INFO": Input columns mapped to key/value fields in the "info" column.

Each entry in "FILTERS" and "FORMAT" sections has the structure

```json
"input column name": {
  "number": 0,       # int
  "type": "str",     # Allowed values: "Flag", "String", "Integer", "Float", 
  "description": "descriptive text"  
}
```

For the "INFO" section, the entries have an additional "new_info_id" field that should contain a name to be used as key in the output "info" column. Note that this string will be repeated in every row of the file, so you may want to keep it short. Mappings from the "INFO" section are only applied if the "new_info_id" value is set.

```json
"input column name": {
  "number": 1,  
  "type": "String",
  "description": "some description",
  "new_info_id": "output attribute name"
}
```

