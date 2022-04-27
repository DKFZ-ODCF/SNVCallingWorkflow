#!/usr/bin/env bash
#
# Copyright (c) 2018 German Cancer Research Center (DKFZ).
#
# Distributed under the MIT License (https://opensource.org/licenses/MIT).
#

set +o

set -ue

module load R/"${RSCRIPT_VERSION:?RSCRIPT_VERSION undefined}"
module load python/"${PYTHON_VERSION:?PYTHON_VERSION undefined}"
module load samtools/"${SAMTOOLS_VERSION:?SAMTOOLS_VERSION undefined}"
module load bcftools/"${BCFTOOLS_VERSION:?BCFTOOLS_VERSION undefined}"
module load htslib/"${HTSLIB_VERSION:?HTSLIB_VERSION undefined}"
module load perl/"${PERL_VERSION:?PERL_VERSION undefined}"
module load bedtools/"${BEDTOOLS_VERSION:?BEDTOOLS_VERSION undefined}"
module load pypy/"${PYPY_VERSION:?PYPY_VERSION undefined}"

echo "Environment is broken for module python/2.7.9. Use pypy or Conda environment" >> /dev/stderr
set +ue
source /dkfz/cluster/virtualenvs/paramasi/python_2.7.9_SNVCalling_pysam_0.16.0.1/bin/activate
set -ue


export BGZIP_BINARY=bgzip
export TABIX_BINARY=tabix
export PERL_BINARY=perl
export PYPY_OR_PYTHON_BINARY="${PYPY_OR_PYTHON_BINARY:-pypy-c}"
export PYTHON_BINARY=python
export RSCRIPT_BINARY=Rscript
export INTERSECTBED_BINARY=intersectBed
export FASTAFROMBED_BINARY=fastaFromBed
export BCFTOOLS_BINARY=bcftools
export BEDTOOLS_BINARY=bedtools
export VCFTOOLS_SORT_BINARY=vcf-sort
export SAMTOOLS_BINARY=samtools
export GHOSTSCRIPT_BINARY=gs
#export ANNOVAR_BINARY=IS SET AS CONFIG VALUE
