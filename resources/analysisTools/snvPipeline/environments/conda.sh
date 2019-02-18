#!/usr/bin/env bash
#
# Copyright (c) 2018 German Cancer Research Center (DKFZ).
#
# Distributed under the MIT License (license terms are at https://github.com/DKFZ-ODCF/IndelCallingWorkflow).
#

# Load a Conda environment.

source activate "${condaEnvironmentName:?No Conda environment name defined. Please set 'condaEnvironmentName'.}" \
    || (echo "Could not load Conda environment '$condaEnvironmentName'" && exit 100)

export BGZIP_BINARY=bgzip
export TABIX_BINARY=tabix
export PERL_BINARY=perl

## This is the slower put cleaner version using CPython/pysam.
## The LSF environment implements a faster version with PyPy/hts-python.
export PYPY_OR_PYTHON_BINARY=python
export PYTHON_BINARY=python
export RSCRIPT_BINARY=Rscript
export INTERSECTBED_BINARY=intersectBed
export BCFTOOLS_BINARY=bcftools
export BEDTOOLS_BINARY=bedtools
export VCFTOOLS_SORT_BINARY=vcf-sort
export SAMTOOLS_BINARY=samtools
export GHOSTSCRIPT_BINARY=gs
export GIT_BINARY=git
#export ANNOVAR_BINARY=IS SET AS CONFIG VALUE

filterPeOverlap() {
    "$PYPY_OR_PYTHON_BINARY" -u "$TOOL_FILTER_PE_OVERLAP" "$@"
}
export -f filterPeOverlap
