#!/bin/bash
#
# Copyright (c) 2018 German Cancer Research Center (DKFZ).
#
# Distributed under the MIT License (https://opensource.org/licenses/MIT).
#

set -xuv

withBamCopy=${withBamCopy-false}
controlBam=${CONTROL_BAMFILE_FULLPATH_BP}
tumorBam=${TUMOR_BAMFILE_FULLPATH_BP}
if [[ $withBamCopy ]]
then
    # Look up lcl scratch first! Consider copying only one file.
    echo "Bam copy not active yet"
fi
export VCF_FOR_SNV_FILES="("
outpath=`dirname ${FILENAME_VCF_RAW}`
for i in ${CHROMOSOME_INDICES[@]}
do
    VCF_FOR_SNV_FILES="${VCF_FOR_SNV_FILES} ${outpath}/snvs_${pid}.${i}.vcf"
done
VCF_FOR_SNV_FILES="${VCF_FOR_SNV_FILES} )"

strlen=`expr ${#VCF_FOR_SNV_FILES} - 2`
allfiles=${VCF_FOR_SNV_FILES:1:$strlen}

# TODO Assumption #C_IND % 2 == 0!
# Take the nth and n - 1
CNT=${#CHROMOSOME_INDICES[@]} # 24 by default
MAX=`expr ${CNT} / 2 - 1` # Arrays are indexed from 0 to (including) 23
for n in `seq 0 $MAX`
do
    export forwardIndex=$n
    export reverseIndex=`expr $CNT - 1 - $n` # Get n - 1

    export CHR_INDEX_FRONT=${CHROMOSOME_INDICES_SORTED[$forwardIndex]}
    export CHR_INDEX_BACK=${CHROMOSOME_INDICES_SORTED[$reverseIndex]}
    export FILE_FRONT=${outpath}/snvs_${pid}.${CHR_INDEX_FRONT}.vcf
    export FILE_BACK=${outpath}/snvs_${pid}.${CHR_INDEX_BACK}.vcf
    export FILE_FRONT_CP=${outpath}/.snvs_${pid}.${CHR_INDEX_FRONT}.vcf_checkpoint
    export FILE_BACK_CP=${outpath}/.snvs_${pid}.${CHR_INDEX_BACK}.vcf_checkpoint
    export FILE_FRONT_CP_OLD=${FILE_FRONT}.cp
    export FILE_BACK_CP_OLD=${FILE_BACK}.cp

    if [[ ! -f ${FILE_FRONT} ]] || [[ ! -f ${FILE_BACK} ]]
    then
    # load snvCalling script twice with a wait after another.
        (FILENAME_VCF_SNVS=${FILE_FRONT} FILENAME_VCF_SNVS_CHECKPOINT=${FILE_FRONT_CP} PARM_CHR_INDEX=${CHR_INDEX_FRONT} FILENAME_VCF_SNVS_CHECKPOINT_OLD=${FILE_FRONT_CP_OLD} bash ${TOOL_SNV_CALLING} ; \
         FILENAME_VCF_SNVS=${FILE_BACK} FILENAME_VCF_SNVS_CHECKPOINT=${FILE_BACK_CP} PARM_CHR_INDEX=${CHR_INDEX_BACK} FILENAME_VCF_SNVS_CHECKPOINT_OLD=${FILE_BACK_CP_OLD} bash ${TOOL_SNV_CALLING}) &
    else
        echo "Both files exist, skipping jobs for ${CHR_INDEX_FRONT}, ${CHR_INDEX_BACK}: (${FILE_FRONT}, ${FILE_BACK}";
    fi
done

wait

# Join!
VCF_FOR_SNV_FILES=${VCF_FOR_SNV_FILES} bash ${TOOL_SNV_JOIN_VCF_FILES} 

