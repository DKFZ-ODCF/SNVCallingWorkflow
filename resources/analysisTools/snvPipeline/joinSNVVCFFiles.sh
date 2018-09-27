#!/bin/bash
#
# Copyright (c) 2018 German Cancer Research Center (DKFZ).
#
# Distributed under the MIT License (https://opensource.org/licenses/MIT).
#
#
#PBS -l walltime=5:00:00
#PBS -l nodes=1:ppn=2
#PBS -l mem=2g
#PBS -m a

set -o pipefail

[[ -z ${TUMOR_BAMFILE_FULLPATH_BP-} ]] && echo "Parameter is missing: TUMOR_BAMFILE_FULLPATH_BP" && exit -10
[[ -z ${FILENAME_VCF_RAW-} ]] && echo "Parameter is missing: FILENAME_VCF_SNVS" && exit -10
[[ -z ${VCF_FOR_SNV_FILES-} ]] && echo "Parameter is missing or empty: VCF_FOR_SNV_FILES" && exit -10

strlen=`expr ${#VCF_FOR_SNV_FILES} - 2`
allfiles=${VCF_FOR_SNV_FILES:1:$strlen}

# Remove snv files after the merge. There will be remaining hidden checkpoint files for the snv calls (for Roddy)
${PERL_BINARY} ${TOOL_HEADERED_FILE_CONCATENATOR} ${allfiles} | ${BGZIP_BINARY} > ${FILENAME_VCF_RAW}.tmp && mv ${FILENAME_VCF_RAW}.tmp ${FILENAME_VCF_RAW} && ${TABIX_BINARY} -p vcf ${FILENAME_VCF_RAW} && [[ ${keepTemporaryVCFFiles-false} == false ]] && rm ${allfiles}
