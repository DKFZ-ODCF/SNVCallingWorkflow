#!/bin/bash
#
# Copyright (c) 2018 German Cancer Research Center (DKFZ).
#
# Distributed under the MIT License (https://opensource.org/licenses/MIT).
#

PLOT_TYPE=${BASE_SCORE_BIAS_PLOT_TYPE:-Differences}

set -xuv
set -o pipefail
[[ -z ${TUMOR_BAMFILE_FULLPATH_BP-} ]] && echo "Parameter is missing: TUMOR_BAMFILE_FULLPATH_BP" && exit -10
[[ -z ${FILENAME_VCF_IN-} ]] && echo "Parameter is missing: FILENAME_VCF_SNVS" && exit -10
[[ -z ${FILENAME_VCF_OUT-} ]] && echo "Parameter is missing: FILENAME_VCF_OUT" && exit -10

declare -r outputDirectory=`dirname ${FILENAME_VCF_OUT}`
relinked=false
# (Re-)link bam file and index. The pipeline cannot handle i.e. .mdup.bam and .mdup.bai
# Look up the index with .bam.bai. If the file exists nothing happens.
if [[ ! -f ${TUMOR_BAMFILE_FULLPATH_BP}.bai ]]
then
    shortBamIdx=`dirname ${TUMOR_BAMFILE_FULLPATH_BP}`"/"`basename ${TUMOR_BAMFILE_FULLPATH_BP} .bam`".bai"
    [[ ! -f ${shortBamIdx} ]] && echo "Bam index file cannot be found!" && exit -15

    # Relink files
    ln -s ${TUMOR_BAMFILE_FULLPATH_BP} ${outputDirectory}
    ln -s ${shortBamIdx} ${outputDirectory}"/"`basename ${TUMOR_BAMFILE_FULLPATH_BP}`".bai"
    # Reset file
    TUMOR_BAMFILE_FULLPATH_BP=${outputDirectory}"/"`basename ${TUMOR_BAMFILE_FULLPATH_BP}`
    relinked=true
fi

scriptDirectory=`dirname ${WRAPPED_SCRIPT}`

# Additional tools
#TOOLS_DIR=${BASEPATH_TOOLS}
FILEBASE=`dirname ${FILENAME_VCF_IN}`/`basename ${FILENAME_VCF_IN} .vcf`

source ${TOOL_ANALYZE_BAM_HEADER}
getRefGenomeAndChrPrefixFromHeader ${TUMOR_BAMFILE_FULLPATH_BP} # Sets CHR_PREFIX and REFERENCE_GENOME


[[ -f ${FILENAME_CHECKPOINT} ]] && rm ${FILENAME_CHECKPOINT}

declare -r filenameSNVVCF=${outputDirectory}/`basename ${FILENAME_VCF_OUT} .gz` # Cut off .gz from the end to get an intermediate file.
declare -r filenameSNVVCFTemp=${filenameSNVVCF}.tmp
declare -r filenameSNVVCFPancan=${outputDirectory}/`basename ${filenameSNVVCF} .vcf`_pancan.vcf.gz
declare -r filenameSNVForAnnovarBed=${FILEBASE}".ForAnnovar.bed"
declare -r filenameSNVForAnnovarBedTemp=${filenameSNVForAnnovarBed}.tmp
declare -r filenameSNVForAnnovarSeqDup=${FILEBASE}".Annovar.seqdup"
declare -r filenameSNVForAnnovarCytoband=${FILEBASE}".Annovar.cytoband"
declare -r filenamePCRerrorMatrix=${outputDirectory}/`basename ${filenameSNVVCF} .vcf`_sequence_error_matrix.txt
declare -r filenameSequencingErrorMatrix=${outputDirectory}/`basename ${filenameSNVVCF} .vcf`_sequencing_error_matrix.txt
declare -r filenameBiasMatrixSeqFile=${outputDirectory}/`basename ${filenameSNVVCF} .vcf`_sequence_specific_bias_matrix.txt
declare -r filenameBiasMatrixSeqingFile=${outputDirectory}/`basename ${filenameSNVVCF} .vcf`_sequencing_specific_bias_matrix.txt
declare -r filenameSomaticSNVsTmp=${outputDirectory}/`basename ${filenameSNVVCF} .vcf`_somatic_snvs_for_bias.vcf
declare -r filenameSequenceErrorPlotPreFilter=${outputDirectory}/`basename ${filenameSNVVCF} .vcf`_sequence_specific_error_plot_before_filter.pdf
declare -r filenameSequencingErrorPlotPreFilter=${outputDirectory}/`basename ${filenameSNVVCF} .vcf`_sequencing_specific_error_plot_before_filter.pdf
declare -r filenameBaseScoreBiasPreFilter=${outputDirectory}/`basename ${filenameSNVVCF} .vcf`_base_score_bias_before_filter.pdf
declare -r filenameSequenceErrorPlotTmp=${outputDirectory}/`basename ${filenameSNVVCF} .vcf`_sequence_specific_error_plot_after_filter_once.pdf
declare -r filenameSequencingErrorPlotTmp=${outputDirectory}/`basename ${filenameSNVVCF} .vcf`_sequencing_specific_error_plot_after_filter_once.pdf
declare -r filenameBaseScoreBiasTmp=${outputDirectory}/`basename ${filenameSNVVCF} .vcf`_base_score_bias_after_filter_once.pdf
declare -r filenameQCValues=${outputDirectory}/`basename ${filenameSNVVCF} .vcf`_QC_values.tsv
outputFilenamePrefix=${outputDirectory}/${SNVFILE_PREFIX}${PID}
declare -r filenameReferenceAlleleBaseScores=${outputFilenamePrefix}_reference_allele_base_qualities.txt.gz
declare -r filenameAlternativeAlleleBaseScores=${outputFilenamePrefix}_alternative_allele_base_qualities.txt.gz
declare -r filenameAlternativeAlleleReadPositions=${outputFilenamePrefix}_alternative_allele_read_positions.txt.gz
declare -r filenameReferenceAlleleReadPositions=${outputFilenamePrefix}_reference_allele_read_positions.txt.gz


declare -r filenamePCRerrorMatrixFirst=${outputDirectory}/`basename ${filenameSNVVCF} .vcf`_sequence_error_matrix_first.txt
declare -r filenameSequencingErrorMatrixFirst=${outputDirectory}/`basename ${filenameSNVVCF} .vcf`_sequencing_error_matrix_first.txt
declare -r filenameBiasMatrixSeqFileFirst=${outputDirectory}/`basename ${filenameSNVVCF} .vcf`_sequence_specific_bias_matrix_first.txt
declare -r filenameBiasMatrixSeqingFileFirst=${outputDirectory}/`basename ${filenameSNVVCF} .vcf`_sequencing_specific_bias_matrix_first.txt

declare -r filenamePCRerrorMatrixSecond=${outputDirectory}/`basename ${filenameSNVVCF} .vcf`_sequence_error_matrix_second.txt
declare -r filenameSequencingErrorMatrixSecond=${outputDirectory}/`basename ${filenameSNVVCF} .vcf`_sequencing_error_matrix_second.txt
declare -r filenameBiasMatrixSeqFileSecond=${outputDirectory}/`basename ${filenameSNVVCF} .vcf`_sequence_specific_bias_matrix_second.txt
declare -r filenameBiasMatrixSeqingFileSecond=${outputDirectory}/`basename ${filenameSNVVCF} .vcf`_sequencing_specific_bias_matrix_second.txt

declare -r filenameControlMedian=${outputDirectory}/`basename ${filenameSNVVCF} .vcf`_control_median.txt

###### Annotate with polymorphisms (dbSNP and 1K genomes) and prepare annovar input file
### Do different things depending on if the file is already zipped or not...

cmdFilter="zcat ${FILENAME_VCF_IN}"

if [[ ${isNoControlWorkflow-false} == "false" ]]; then
    cmdFilter="${cmdFilter} | perl ${TOOL_MEDIAN} - ${filenameControlMedian}"
fi

cmdFilter="${cmdFilter} | ${PERL_BINARY} ${TOOL_ANNOTATE_VCF_FILE} --tabix_bin=${TABIX_BINARY} -a - -b ${DBSNP} --columnName=${DBSNP_COL} --reportMatchType  --bAdditionalColumn=2 --reportLevel 4 | \
    ${PERL_BINARY} ${TOOL_ANNOTATE_VCF_FILE} --tabix_bin=${TABIX_BINARY} -a - -b ${KGENOME} --columnName=${KGENOMES_COL}  --reportMatchType --bAdditionalColumn=2 --reportLevel 4 | \
    ${PERL_BINARY} ${TOOL_ANNOTATE_VCF_FILE} --tabix_bin=${TABIX_BINARY} -a - -b ${ExAC} --columnName ${ExAC_COL} --bFileType vcf --reportMatchType --reportLevel 4 | \
    ${PERL_BINARY} ${TOOL_ANNOTATE_VCF_FILE} --tabix_bin=${TABIX_BINARY} -a - -b ${EVS} --columnName ${EVS_COL} --bFileType vcf --reportMatchType --reportLevel 4 | \
    ${PERL_BINARY} ${TOOL_ANNOTATE_VCF_FILE} -a - -b ${GNOMAD_WES_ALL_SNV} --columnName=${GNOMAD_WES_COL} --bFileType vcf --reportMatchType --reportLevel 4 | \
    ${PERL_BINARY} ${TOOL_ANNOTATE_VCF_FILE} -a - -b ${GNOMAD_WGS_ALL_SNV} --columnName=${GNOMAD_WGS_COL} --bFileType vcf --reportMatchType --reportLevel 4 | \
    ${PERL_BINARY} ${TOOL_ANNOTATE_VCF_FILE} --tabix_bin=${TABIX_BINARY} -a - -b ${LOCALCONTROL_WGS} --columnName ${LOCALCONTROL_WGS_COL} --bFileType vcf --reportMatchType --reportLevel 4 | \
    ${PERL_BINARY} ${TOOL_ANNOTATE_VCF_FILE} --tabix_bin=${TABIX_BINARY} -a - -b ${LOCALCONTROL_WES} --columnName ${LOCALCONTROL_WES_COL} --bFileType vcf --reportMatchType --reportLevel 4"

if [[ -f ${RECURRENCE} ]]; then
    cmdFilter="${cmdFilter} | ${PERL_BINARY} ${TOOL_ANNOTATE_VCF_FILE} -a - -b ${RECURRENCE} --columnName ${RECURRENCE_COL} --bFileType vcf"
fi
cmdFilter="${cmdFilter} | tee ${filenameSNVVCFTemp} | perl ${TOOL_VCF_TO_ANNOVAR} ${CHR_PREFIX} ${CHR_SUFFIX} > ${filenameSNVForAnnovarBedTemp}"

eval ${cmdFilter}

[[ "$?" != 0 ]] && echo "There was a non-zero exit code in the polymorphism annotation pipe; exiting..." && exit 2

mv ${filenameSNVVCFTemp} ${filenameSNVVCF}
mv ${filenameSNVForAnnovarBedTemp} ${filenameSNVForAnnovarBed}


###### Basic annotation with Annovar
# Gene annotation with annovar
ANNOVAR_FILE=${ANNOVAR_BUILDVER}_`eval echo $ANNOVAR_DBTYPE | cut -d " " -f 2`.txt
ANNOVAR_DBFILEPATH=${ANNOVAR_DBPATH}/${ANNOVAR_FILE}
[[ ! -f ${ANNOVAR_DBFILEPATH} ]]  && echo "Gene annotation database not found. Check ANNOVAR_DBTYPE." && exit -16
[[ ${runGeneAnnovar-true} == "true" ]] && ${ANNOVAR_BINARY} --buildver=${ANNOVAR_BUILDVER} ${ANNOVAR_DBTYPE} ${filenameSNVForAnnovarBed} ${ANNOVAR_DBPATH}

# segdup annotation with annovar
${ANNOVAR_BINARY} --buildver=${ANNOVAR_BUILDVER} -regionanno -dbtype segdup --outfile=${filenameSNVForAnnovarSeqDup} ${filenameSNVForAnnovarBed} ${ANNOVAR_DBPATH}
av_segdup=`ls ${filenameSNVForAnnovarSeqDup}*genomicSuperDups`

# cytoband annotation with annovar
if [[ ${runCytoband-true} == "true" ]]
then
	${ANNOVAR_BINARY} --buildver=$ANNOVAR_BUILDVER -regionanno -dbtype band --outfile=${filenameSNVForAnnovarCytoband} ${filenameSNVForAnnovarBed} ${ANNOVAR_DBPATH}
	av_cytoband=`ls ${filenameSNVForAnnovarCytoband}*cytoBand`
fi


if [[ ${runGeneAnnovar-true} == "true" && ${runCytoband-true} == "true" ]]
then
    perl ${TOOL_CONVERT_NEWCOLS_TO_VCF} --vcfFile=${filenameSNVVCF} \
        --newColFile="${PERL_BINARY} ${TOOL_PROCESS_ANNOVAR} ${filenameSNVForAnnovarBed}.variant_function ${filenameSNVForAnnovarBed}.exonic_variant_function |" \
        --newColHeader=$ANNOVAR_GENEANNO_COLS --chrPrefix=$CHR_PREFIX --chrSuffix=$CHR_SUFFIX --reportColumns="3,4,5,6" --bChrPosEnd="0,1,2" |
    perl ${TOOL_CONVERT_NEWCOLS_TO_VCF} --vcfFile="-" --newColFile=${av_segdup} --newColHeader=$ANNOVAR_SEGDUP_COL --chrPrefix=$CHR_PREFIX --chrSuffix=$CHR_SUFFIX \
        --reportColumns="1" --bChrPosEnd="2,7,8"  |
    perl ${TOOL_CONVERT_NEWCOLS_TO_VCF} --vcfFile="-" --newColFile=${av_cytoband} --newColHeader=$ANNOVAR_CYTOBAND_COL --chrPrefix=$CHR_PREFIX --chrSuffix=$CHR_SUFFIX \
        --reportColumns="1" --bChrPosEnd="2,7,8" > ${filenameSNVVCFTemp}
elif [[ ${runGeneAnnovar-true} == "true" ]]
then
    perl ${TOOL_CONVERT_NEWCOLS_TO_VCF} --vcfFile=${filenameSNVVCF} \
        --newColFile="${PERL_BINARY} ${TOOL_PROCESS_ANNOVAR} ${filenameSNVForAnnovarBed}.variant_function ${filenameSNVForAnnovarBed}.exonic_variant_function |" \
        --newColHeader=$ANNOVAR_GENEANNO_COLS --chrPrefix=$CHR_PREFIX --chrSuffix=$CHR_SUFFIX --reportColumns="3,4,5,6" --bChrPosEnd="0,1,2" |
    perl ${TOOL_CONVERT_NEWCOLS_TO_VCF} --vcfFile="-" --newColFile=${av_segdup} --newColHeader=$ANNOVAR_SEGDUP_COL --chrPrefix=$CHR_PREFIX --chrSuffix=$CHR_SUFFIX \
        --reportColumns="1" --bChrPosEnd="2,7,8" > ${filenameSNVVCFTemp}
else
    perl ${TOOL_CONVERT_NEWCOLS_TO_VCF} --vcfFile="${filenameSNVVCF}" --newColFile=${av_segdup} --newColHeader=$ANNOVAR_SEGDUP_COL --chrPrefix=$CHR_PREFIX --chrSuffix=$CHR_SUFFIX \
        --reportColumns="1" --bChrPosEnd="2,7,8" > ${filenameSNVVCFTemp}
fi
exitCode=$?
[[ $exitCode == 0 ]] && mv ${filenameSNVVCFTemp} ${filenameSNVVCF}
[[ $exitCode != 0 ]] && echo "There was a non-zero exit code in the Annovar annotation pipe; temp file ${filenameSNVVCFTemp} not moved back" && exit 2

snv_reliability_pipe=`perl ${TOOL_CREATEPIPES} ${filenameSNVVCF} ${PARAMETER_FILE} ${TOOL_ANNOTATE_VCF_FILE} SNV_RELIABILITY ${TABIX_BINARY}`

if [[ "$?" != 0 ]] || [[ -z "$snv_reliability_pipe" ]]; then echo "problem when generating SNV_RELIABILITY pipe. Exiting..."; exit 2; fi

eval $snv_reliability_pipe > ${filenameSNVVCFTemp}
exitCode=$?
[[ $exitCode == 0 ]] && mv ${filenameSNVVCFTemp} ${filenameSNVVCF}
[[ $exitCode != 0 ]] && echo "There was a non-zero exit code in the SNV_RELIABILITY pipe; temp file ${filenameSNVVCFTemp} not moved back" && exit 2

# confidence annotation with/without germline
if [[ ${isNoControlWorkflow-false} == "false" ]]
then
    MAX_CONTROL_COV=0
    MAX_CONTROL_COV=`cat ${filenameControlMedian}`
    [[ ${MAX_CONTROL_COV} -lt 1 ]] && echo "Median was not calculated correctly" && exit 3
    rm ${filenameControlMedian}
    [[ ${CONFIDENCE_OPTS} != *"-t "* ]] && [[ ${CONFIDENCE_OPTS} != *"-t= "* ]] && [[ ${CONFIDENCE_OPTS} != *"--cutoff "* ]] && [[ ${CONFIDENCE_OPTS} != *"--cutoff= "* ]] && CONFIDENCE_OPTS=${CONFIDENCE_OPTS}" -t ${MAX_CONTROL_COV}"
    noControlFlag=""
else
    noControlFlag="--nocontrol"
fi

# for somatic SNVs only: if there are overlapping reads (same name) at the position, count the base only once;
# also remove counts for bases that are neither reference nor variant, and change DP4 field accordingly
# this needs the BAM file to perform a pileup

eval "tumorbamfullpath=$TUMOR_BAMFILE_FULLPATH_BP"           # $RESULTS_PER_PIDS_DIR/${pid}/$ALIGNMENT_SUBDIR/$tumorbam
# and the basequal and mapqual used for mpileup-bcftools.
# if not specified, set to defaults 30 and 13, respectively
mapqual=`echo $MPILEUP_OPTS | perl -ne '($qual) = $_ =~ /\-q\s*(\d+)/;print $qual'`
mapqual=${mapqual:-30}
basequal=`echo $MPILEUP_OPTS | perl -ne '($qual) = $_ =~ /\-Q\s*(\d+)/;print $qual'`
basequal=${basequal:-13}

if [[ ${runOnPancan-false} == true ]]; then
    CONFIDENCE_OPTS_PANCAN=${CONFIDENCE_OPTS}" -o ${filenameSNVVCFPancan}"
else
    CONFIDENCE_OPTS_PANCAN=${CONFIDENCE_OPTS}
fi

npConfidence=${RODDY_SCRATCH}/snvAnnotationFIFO
mkfifo ${npConfidence}
if [[ ${isNoControlWorkflow-false} == "false" ]]; then
    cat ${filenameSNVVCF} > ${npConfidence} &
else
    cat ${filenameSNVVCF} | ${PYPY_OR_PYTHON_BINARY} -u ${TOOL_CONFIDENCE_ANNOTATION} ${noControlFlag} -i - ${CONFIDENCE_OPTS} -a 0 > ${npConfidence} &
fi

# create BaseScore FIFOs and their consumer processes (zip and write to target file)
# BaseScore FIFOS will be filled by ${TOOL_FILTER_PE_OVERLAP}
mkfifo ${filenameAlternativeAlleleBaseScores}_NP ${filenameReferenceAlleleBaseScores}_NP ${filenameAlternativeAlleleReadPositions}_NP ${filenameReferenceAlleleReadPositions}_NP
cat ${filenameAlternativeAlleleBaseScores}_NP | ${BGZIP_BINARY} -f >${filenameAlternativeAlleleBaseScores} & zipAlternativeAlleleBaseScores=$!
cat ${filenameReferenceAlleleBaseScores}_NP | ${BGZIP_BINARY} -f >${filenameReferenceAlleleBaseScores} & zipReferenceAlleleBaseScores=$!
cat ${filenameAlternativeAlleleReadPositions}_NP | ${BGZIP_BINARY} -f >${filenameAlternativeAlleleReadPositions} & zipAlternativeAlleleReadPositions=$!
cat ${filenameReferenceAlleleReadPositions}_NP | ${BGZIP_BINARY} -f >${filenameReferenceAlleleReadPositions} & zipReferenceAlleleReadPositions=$!

# If this is for the pancancer workflow, then also create a DKFZ specific file.
if [[ ${runArtifactFilter-true} == true ]]
then
	cat ${npConfidence} | filterPeOverlap ${noControlFlag} --alignmentFile=${tumorbamfullpath} --mapq=$mapqual --baseq=$basequal --qualityScore=phred --maxNumberOfMismatchesInRead=${NUMBER_OF_MISMACTHES_THRESHOLD--1} --altBaseQualFile=${filenameAlternativeAlleleBaseScores}_NP --refBaseQualFile=${filenameReferenceAlleleBaseScores}_NP --altBasePositionsFile=${filenameAlternativeAlleleReadPositions}_NP --refBasePositionsFile ${filenameReferenceAlleleReadPositions}_NP --referenceFile=${REFERENCE_GENOME} \
						| ${PYPY_OR_PYTHON_BINARY} -u ${TOOL_CONFIDENCE_ANNOTATION} ${noControlFlag} -i - ${CONFIDENCE_OPTS} -a 0 -f ${filenameSomaticSNVsTmp} > ${filenameSNVVCFTemp}.tmp

    [[ $? != 0 ]] && echo "Error in first iteration of confidence annotation" && exit 2

	NRSOMSNV=`grep -v "^#" ${filenameSomaticSNVsTmp} | wc -l`
	echo -e "SOMATIC_SNVS_UNFILTERED\t${NRSOMSNV}">> ${filenameQCValues}

	mv ${filenameSNVVCFTemp}.tmp ${filenameSNVVCFTemp}

    wait ${zipAlternativeAlleleBaseScores} ; [[ $? -gt 0 ]] && echo "Error from zipAlternativeAlleleBaseScores" && exit 31
    wait ${zipReferenceAlleleBaseScores} ; [[ $? -gt 0 ]] && echo "Error from zipReferenceAlleleBaseScores" && exit 32
    wait ${zipAlternativeAlleleReadPositions} ; [[ $? -gt 0 ]] && echo "Error from zipAlternativeAlleleReadPositions" && exit 33
    wait ${zipReferenceAlleleReadPositions} ; [[ $? -gt 0 ]] && echo "Error from zipReferenceAlleleReadPositions" && exit 34
    rm ${filenameAlternativeAlleleBaseScores}_NP ${filenameReferenceAlleleBaseScores}_NP ${filenameAlternativeAlleleReadPositions}_NP ${filenameReferenceAlleleReadPositions}_NP

	${PYTHON_BINARY} ${TOOL_CREATE_ERROR_PLOTS} --vcfFile=${filenameSomaticSNVsTmp} --referenceFile=NA --outputFile=${filenameSequencingErrorPlotPreFilter} --errorType=sequencing_specific --errorFile=${filenameSequencingErrorMatrix} --plot_title='Sequencing strand bias before guanine oxidation filter'

	[[ $? != 0 ]] && echo "Error in first creation of error matrix and plot (sequencing)" && exit 3

	${PYTHON_BINARY} ${TOOL_CREATE_ERROR_PLOTS} --vcfFile=${filenameSomaticSNVsTmp} --referenceFile=NA --outputFile=${filenameSequenceErrorPlotPreFilter} --errorType=sequence_specific --errorFile=${filenamePCRerrorMatrix} --plot_title='PCR strand bias before guanine oxidation filter'

	[[ $? != 0 ]] && echo "Error in first creation of error matrix and plot (sequence/PCR)" && exit 4

    if [[ ${generateExtendedQcPlots} == true ]]; then
        cp ${filenameSomaticSNVsTmp} ${filenameSomaticSNVsTmp}.forBSBPlotsBeforeFiltering
        ${RSCRIPT_BINARY} ${TOOL_PLOT_BASE_SCORE_BIAS} -v ${filenameSomaticSNVsTmp}.forBSBPlotsBeforeFiltering -r ${filenameReferenceAlleleBaseScores} -a ${filenameAlternativeAlleleBaseScores} -t ${basequal} -p ${PLOT_TYPE} -o ${filenameBaseScoreBiasPreFilter} -d "Base Quality Bias Plot for PID ${PID} before guanine oxidation filter" & plotBaseScoreBiasBeforeFiltering=$!
    fi

	cat ${filenameSNVVCFTemp} | ${PYTHON_BINARY} ${TOOL_FLAG_BIAS} --vcfFile="/dev/stdin" --referenceFile=${REFERENCE_GENOME} --sequence_specificFile=${filenamePCRerrorMatrix} --sequencing_specificFile=${filenameSequencingErrorMatrix} --numReads=${nReads} --numMuts=${nMuts} --biasPValThreshold=${biasPValThreshold} --biasRatioThreshold=${biasRatioThreshold} --biasRatioMinimum=${biasRatioMinimum} --maxNumOppositeReadsSequencingWeakBias=${maxNumOppositeReadsSequencingWeakBias} --maxNumOppositeReadsSequenceWeakBias=${maxNumOppositeReadsSequenceWeakBias} --maxNumOppositeReadsSequencingStrongBias=${maxNumOppositeReadsSequencingStrongBias} --maxNumOppositeReadsSequenceStrongBias=${maxNumOppositeReadsSequenceStrongBias} --ratioVcf=${rVcf} --bias_matrixSeqFile=${filenameBiasMatrixSeqFile} --bias_matrixSeqingFile=${filenameBiasMatrixSeqingFile} --vcfFileFlagged="/dev/stdout" | \
	${PYPY_OR_PYTHON_BINARY} -u ${TOOL_CONFIDENCE_ANNOTATION} ${noControlFlag} -i - ${CONFIDENCE_OPTS} -a 1 -f ${filenameSomaticSNVsTmp} > ${filenameSNVVCFTemp}.tmp

	[[ $? != 0 ]] && echo "Error in first filtering and/or second interation of confidence annotation" && exit 5

	mv ${filenameSNVVCFTemp}.tmp ${filenameSNVVCFTemp}
	mv ${filenamePCRerrorMatrix} ${filenamePCRerrorMatrixFirst}
	mv ${filenameSequencingErrorMatrix} ${filenameSequencingErrorMatrixFirst}
	mv ${filenameBiasMatrixSeqFile} ${filenameBiasMatrixSeqFileFirst}
	mv ${filenameBiasMatrixSeqingFile} ${filenameBiasMatrixSeqingFileFirst}

	${PYTHON_BINARY} ${TOOL_CREATE_ERROR_PLOTS} --vcfFile=${filenameSomaticSNVsTmp} --referenceFile=NA --outputFile=${filenameSequencingErrorPlotTmp} --errorType=sequencing_specific --errorFile=${filenameSequencingErrorMatrix} --plot_title='Sequencing strand bias after first round of guanine oxidation filter'

    [[ $? != 0 ]] && echo "Error in second creation of error matrix and plot (sequencing)" && exit 6

	${PYTHON_BINARY} ${TOOL_CREATE_ERROR_PLOTS} --vcfFile=${filenameSomaticSNVsTmp} --referenceFile=NA --outputFile=${filenameSequenceErrorPlotTmp} --errorType=sequence_specific --errorFile=${filenamePCRerrorMatrix} --plot_title='PCR strand bias after first round of guanine oxidation filter'

	[[ $? != 0 ]] && echo "Error in second creation of error matrix and plot (sequence/PCR)" && exit 7

    if [[ ${generateExtendedQcPlots} == true ]]; then
        cp ${filenameSomaticSNVsTmp} ${filenameSomaticSNVsTmp}.forBSBPlotsAfter1stFiltering
        ${RSCRIPT_BINARY} ${TOOL_PLOT_BASE_SCORE_BIAS} -v ${filenameSomaticSNVsTmp}.forBSBPlotsAfter1stFiltering -r ${filenameReferenceAlleleBaseScores} -a ${filenameAlternativeAlleleBaseScores} -t ${basequal} -p ${PLOT_TYPE} -o ${filenameBaseScoreBiasTmp} -d "Base Quality Bias Plot for PID ${PID} after first round of guanine oxidation filter" & plotBaseScoreBiasAfterFirstFiltering=$!
    fi

	cat ${filenameSNVVCFTemp} | ${PYTHON_BINARY} -u ${TOOL_FLAG_BIAS} --vcfFile="/dev/stdin" --referenceFile=${REFERENCE_GENOME} --sequence_specificFile=${filenamePCRerrorMatrix} --sequencing_specificFile=${filenameSequencingErrorMatrix} --numReads=${nReads} --numMuts=${nMuts} --biasPValThreshold=${biasPValThreshold} --biasRatioThreshold=${biasRatioThreshold} --biasRatioMinimum=${biasRatioMinimum} --maxNumOppositeReadsSequencingWeakBias=${maxNumOppositeReadsSequencingWeakBias} --maxNumOppositeReadsSequenceWeakBias=${maxNumOppositeReadsSequenceWeakBias} --maxNumOppositeReadsSequencingStrongBias=${maxNumOppositeReadsSequencingStrongBias} --maxNumOppositeReadsSequenceStrongBias=${maxNumOppositeReadsSequenceStrongBias} --ratioVcf=${rVcf} --bias_matrixSeqFile=${filenameBiasMatrixSeqFile} --bias_matrixSeqingFile=${filenameBiasMatrixSeqingFile} --vcfFileFlagged="/dev/stdout" | \
	${PYPY_OR_PYTHON_BINARY} -u ${TOOL_CONFIDENCE_ANNOTATION} ${noControlFlag} -i - ${CONFIDENCE_OPTS_PANCAN} -a 2  > ${filenameSNVVCF}

	[[ $? != 0 ]] && echo "Error in second filtering and/or third iteration of confidence annotation" && exit 8

	mv ${filenamePCRerrorMatrix} ${filenamePCRerrorMatrixSecond}
	mv ${filenameSequencingErrorMatrix} ${filenameSequencingErrorMatrixSecond}
	mv ${filenameBiasMatrixSeqFile} ${filenameBiasMatrixSeqFileSecond}
	mv ${filenameBiasMatrixSeqingFile} ${filenameBiasMatrixSeqingFileSecond}

	rm ${filenameSomaticSNVsTmp}
	rm ${filenameSNVVCFTemp}

	[[ $? != 0 ]] && echo "Error in moving the vcf file and index or in removing the temporary files" && exit 9

	if [[ ${generateExtendedQcPlots} == true ]]; then
    wait ${plotBaseScoreBiasBeforeFiltering} ; [[ $? -gt 0 ]] && echo "Error in first creation of base score bias plot" && exit 37
    [[ -f ${filenameSomaticSNVsTmp}.forBSBPlotsBeforeFiltering ]] && rm ${filenameSomaticSNVsTmp}.forBSBPlotsBeforeFiltering
    wait ${plotBaseScoreBiasAfterFirstFiltering} ; [[ $? -gt 0 ]] && echo "Error in second creation of base score bias plot" && exit 38
    [[ -f ${filenameSomaticSNVsTmp}.forBSBPlotsAfter1stFiltering ]] && rm ${filenameSomaticSNVsTmp}.forBSBPlotsAfter1stFiltering
	fi
else
	cat ${npConfidence} | filterPeOverlap ${noControlFlag} --alignmentFile=${tumorbamfullpath} --mapq=$mapqual --baseq=$basequal --qualityScore=phred --maxNumberOfMismatchesInRead=${NUMBER_OF_MISMACTHES_THRESHOLD--1} --altBaseQualFile=${filenameAlternativeAlleleBaseScores}_NP --refBaseQualFile=${filenameReferenceAlleleBaseScores}_NP --altBasePositionsFile=${filenameAlternativeAlleleReadPositions}_NP --refBasePositionsFile=${filenameReferenceAlleleReadPositions}_NP --referenceFile=${REFERENCE_GENOME} | \
	${PYPY_OR_PYTHON_BINARY} -u ${TOOL_CONFIDENCE_ANNOTATION} ${noControlFlag} -i - ${CONFIDENCE_OPTS_PANCAN} > ${filenameSNVVCFTemp}

    exitCode=$?
    [[ $exitCode == 0 ]] && [[ -f ${filenameSNVVCFTemp} ]] && mv ${filenameSNVVCFTemp} ${filenameSNVVCF}
    [[ $exitCode != 0 ]] && echo "SNV confidenceAnnotation with germline pipe returned non-zero exit code; temp file ${filenameSNVVCFTemp} not moved back" && exit 21
    #${BGZIP_BINARY} -f ${filenameSNVVCF} && ${TABIX_BINARY} -f -p vcf ${FILENAME_VCF_OUT}

    wait ${zipAlternativeAlleleBaseScores} ; [[ $? -gt 0 ]] && echo "Error from zipAlternativeAlleleBaseScores" && exit 31
    wait ${zipReferenceAlleleBaseScores} ; [[ $? -gt 0 ]] && echo "Error from zipReferenceAlleleBaseScores" && exit 32
    wait ${zipAlternativeAlleleReadPositions} ; [[ $? -gt 0 ]] && echo "Error from zipAlternativeAlleleReadPositions" && exit 33
    wait ${zipReferenceAlleleReadPositions} ; [[ $? -gt 0 ]] && echo "Error from zipReferenceAlleleReadPositions" && exit 34

    rm ${filenameAlternativeAlleleBaseScores}_NP ${filenameReferenceAlleleBaseScores}_NP ${filenameAlternativeAlleleReadPositions}_NP ${filenameReferenceAlleleReadPositions}_NP
fi
rm ${npConfidence}

#cp ${filenameSNVVCF} ${filenameSNVVCFTemp}.afterFilter2.vcf
${BGZIP_BINARY} -f ${filenameSNVVCF} && ${TABIX_BINARY} -f -p vcf ${FILENAME_VCF_OUT}
[[ $? != 0 ]] && echo "Error in creation of bgzipped vcf file and tabix index for it" && exit 41

# If this is for the pancancer workflow, then also zip away the DKFZ only file.
[[ -f ${filenameSNVVCFPancan} ]] && ${TABIX_BINARY} -f -p vcf ${filenameSNVVCFPancan}

annofile=$(ls ${FILEBASE}*Annovar* 2> /dev/null | wc -l)
[[ "$annofile" != 0 ]] && rm ${FILEBASE}*Annovar*

[[ ${relinked} == true ]] && rm ${TUMOR_BAMFILE_FULLPATH_BP}.bai && rm ${TUMOR_BAMFILE_FULLPATH_BP}

touch ${FILENAME_CHECKPOINT}
