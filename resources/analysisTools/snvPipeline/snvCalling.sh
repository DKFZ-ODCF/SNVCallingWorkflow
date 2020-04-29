#!/bin/bash
#
# Copyright (c) 2018 German Cancer Research Center (DKFZ).
#
# Distributed under the MIT License (https://opensource.org/licenses/MIT).
#
#PBS -l walltime=20:00:00,nodes=1

set -o pipefail

# Create variables and check input parameters
CHR_PREFIX=${CHR_PREFIX-}
CHR_SUFFIX=${CHR_SUFFIX-}
TUMOR_BAMFILE_FULLPATH_BP=${TUMOR_BAMFILE_FULLPATH_BP-}
CONTROL_BAMFILE_FULLPATH_BP=${CONTROL_BAMFILE_FULLPATH_BP-}
FILENAME_VCF_SNVS=${FILENAME_VCF_SNVS-}

[[ -z $TUMOR_BAMFILE_FULLPATH_BP ]] && echo "Parameter is missing: TUMOR_BAMFILE_FULLPATH_BP" && exit -10
if [[ ${isNoControlWorkflow-false} == "false" ]] ; then
    [[ -z $CONTROL_BAMFILE_FULLPATH_BP ]] && echo "Parameter is missing: CONTROL_BAMFILE_FULLPATH_BP" && exit -10
fi
[[ -z $FILENAME_VCF_SNVS ]] && echo "Parameter is missing: FILENAME_VCF_SNVS" && exit -10

PARM_CHR_INDEX=${PARM_CHR_INDEX-}

# Ignore chr prefix and set it manually
source ${TOOL_ANALYZE_BAM_HEADER}
getRefGenomeAndChrPrefixFromHeader ${TUMOR_BAMFILE_FULLPATH_BP} # Sets CHR_PREFIX and REFERENCE_GENOME

chr=""
if [[ -n "$PARM_CHR_INDEX" ]]
then
    chr=${CHR_PREFIX}${PARM_CHR_INDEX}${CHR_SUFFIX}
else
    chr=${CHR_PREFIX}${CHROMOSOME_INDICES[$((RODDY_JOBARRAYINDEX-1))]}${CHR_SUFFIX}
    # File names can contain X and Y. The pattern contains RODDY_JOBARRAYINDEX so this has to be changed temporary as the value is always numeric!
#    temp=$RODDY_JOBARRAYINDEX
#    RODDY_JOBARRAYINDEX=${CHROMOSOME_INDICES[$((temp-1))]}
fi
#env > ${FILENAME_VCF_SNVS}.${chr}.env

runCompareGermline=${runCompareGermline-true} # Set to true if it is not set


# Substitute variables in passed parameters
#FILENAME_VCF_SNVS=`eval echo ${FILENAME_VCF_SNVS//#/$}`

#[[ -n "$PARM_CHR_INDEX" ]] && RODDY_JOBARRAYINDEX=$temp

scriptDirectory=`dirname ${WRAPPED_SCRIPT}`

# Additional tools
TOOLS_DIR=`dirname ${TOOL_SEQ_CONTEXT_ANNOTATOR}`


tumorbamfullpath=$TUMOR_BAMFILE_FULLPATH_BP
tumor_basename=`basename ${tumorbamfullpath} .bam`

MPILEUP_SUBDIR=`dirname ${FILENAME_VCF_SNVS}`
SCRATCH_DIR=${RODDY_SCRATCH}

filenameMPileupTemp="${SCRATCH_DIR}/$tumor_basename.${chr}.tmp"
filenameMPileupTempError="${MPILEUP_SUBDIR}/${MPILEUPOUT_PREFIX}${PID}.${chr}.mpileup_error.tmp"
filenameMPileupOut="${SCRATCH_DIR}/${MPILEUPOUT_PREFIX}${PID}.${chr}.vcf"
filenameMPileupError="${MPILEUP_SUBDIR}/${MPILEUPOUT_PREFIX}${PID}.${chr}.error"

FILENAME_VCF_SNVS_TEMP=${FILENAME_VCF_SNVS}.tmp

#[[ ! -f $filenameMPileupTemp ]] && 
${BCFTOOLS_BINARY} mpileup ${MPILEUP_OPTS} -r ${chr} -f ${REFERENCE_GENOME} $tumorbamfullpath | ${BCFTOOLS_BINARY} call ${BCFTOOLS_OPTS} > $filenameMPileupTemp

if [[ "$?" != 0 ]]
then
	echo "There was a non-zero exit code in the mpileup / bcftools pipe; copying mpileup tmp to snv subdir and exiting..."
	mv ${filenameMPileupTemp} ${filenameMPileupTempError}
	exit 2
fi

if [[ ${isNoControlWorkflow-false} == "false" ]]; then
    snvOut=${filenameMPileupOut}
else
    snvOut=${FILENAME_VCF_SNVS_TEMP}
fi

firstLineVCF=`cat $filenameMPileupTemp | grep -v "^#" | head -n1`
if [[ -n "$firstLineVCF" ]]; then
    # filter out mutations with strand bias next to the motif GGnnGG (or CCnnCC if snv is on reverse strand)
    ${PERL_BINARY} "$TOOL_SEQ_CONTEXT_ANNOTATOR" "$FASTAFROMBED_BINARY" "$filenameMPileupTemp" "$REFERENCE_GENOME" 10 "$TOOLS_DIR" | \
        ${PYTHON_BINARY} "$TOOL_RAW_SNV_FILTER" --outf="$snvOut" "$RAW_SNV_FILTER_OPTIONS" # --inf=$MPILEUP_SUBDIR/forStrandBiasFilter.${chr}.bed

    if [[ "$?" == 0 ]]
    then
	    rm ${filenameMPileupTemp}
    else
	    echo "There was a non-zero exit code in the seqContext_annotator / rawSnvFilter pipe; copying mpileup tmp to snv subdir and exiting..."
	    mv ${filenameMPileupTemp} ${filenameMPileupTempError}
	    exit 2
    fi

    if [[ ${isNoControlWorkflow-false} == "false" ]]; then
        if [[ ${runCompareGermline} == true ]]; then
            NP_MPILEUP=${SCRATCH_DIR}/NP_MPILEUP_CHR${chr}
            mkfifo $NP_MPILEUP

            # if there is a germline BAM, first look up these positions in the control file by just piling up the bases, this is NO SNV calling
            ${SAMTOOLS_BINARY} mpileup ${MPILEUPCONTROL_OPTS} -r ${chr} -l ${filenameMPileupOut} -f ${REFERENCE_GENOME} ${CONTROL_BAMFILE_FULLPATH_BP} > ${NP_MPILEUP} &
            ${PERL_BINARY} ${TOOL_VCF_PILEUP_COMPARE} ${filenameMPileupOut} $NP_MPILEUP "Header" > ${FILENAME_VCF_SNVS_TEMP}

            rm $NP_MPILEUP

            if [[ "$?" == 0  ]]; then
                rm ${filenameMPileupOut}
            else
                echo "vcf_pileup_compare pipe returned non-zero exit code; not moving tmp output files  (${FILENAME_VCF_SNVS_TEMP} to final name; moving $filenameMPileupOut to $MPILEUP_SUBDIR/"
                mv ${filenameMPileupOut} ${filenameMPileupError}
                exit 2
            fi
        fi
    fi
    mv ${FILENAME_VCF_SNVS_TEMP} ${FILENAME_VCF_SNVS}
else
    echo "No SNV was found in ${chr}."
    awk '{print $0"\tSEQUENCE_CONTEXT"}' ${filenameMPileupTemp} > ${FILENAME_VCF_SNVS}
    rm ${filenameMPileupTemp}
fi

touch ${FILENAME_VCF_SNVS_CHECKPOINT}
[[ ${FILENAME_VCF_SNVS_CHECKPOINT_OLD} ]] && touch ${FILENAME_VCF_SNVS_CHECKPOINT_OLD}

# Cleanup
# [[ ${useCustomScratchDir} == true ]] && rm -rf ${SCRATCH_DIR}

exit 0
