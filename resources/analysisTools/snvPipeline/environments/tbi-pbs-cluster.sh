#!/usr/bin/env bash

module load R/"${RSCRIPT_VERSION:?RSCRIPT_VERSION undefined}"
module load python/"${PYTHON_VERSION:?PYTHON_VERSION undefined}"
module load samtools/"${SAMTOOLS_VERSION:?SAMTOOLS_VERSION undefined}"
module load htslib/"${HTSLIB_VERSION:?HTSLIB_VERSION undefined}"
module load perl/"${PERL_VERSION:?PERL_VERSION undefined}"
module load bedtools/"${BEDTOOLS_VERSION:?BEDTOOLS_VERSION undefined}"

source /ibios/tbi_cluster/virtualenvs/warsow/python_2.7.9_SNVCalling_1.2.166-1_PBS/bin/activate

export BGZIP_BINARY=bgzip
export TABIX_BINARY=tabix
export INTERSECTBED_BINARY=intersectBed
export FASTAFROMBED_BINARY=fastaFromBed
export BCFTOOLS_BINARY=bcftools
export BEDTOOLS_BINARY=bedtools
export VCFTOOLS_SORT_BINARY=vcf-sort
export SAMTOOLS_BINARY=samtools
#export ANNOVAR_BINARY=IS SET AS CONFIG VALUE

export PERL_BINARY=perl
export PYTHON_BINARY=python
export RSCRIPT_BINARY=Rscript

export GHOSTSCRIPT_BINARY=gs

