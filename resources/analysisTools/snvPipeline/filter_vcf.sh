#!/bin/bash

source ${CONFIG_FILE}

PLOT_TYPE=${BASE_SCORE_BIAS_PLOT_TYPE:-Differences}

set -o pipefail

[[ -f ${FILENAME_CHECKPOINT} ]] && rm ${FILENAME_CHECKPOINT}
mpileupDirectory=`dirname ${FILENAME_VCF}`

RUN_PLOTS=${RUN_PLOTS-1}
RUN_PUREST=${RUN_PUREST-1}
CONVERT_BINARY=${CONVERT_BINARY-convert}

source ${TOOL_ANALYZE_BAM_HEADER}
getRefGenomeAndChrPrefixFromHeader ${TUMOR_BAMFILE_FULLPATH_BP} # Sets CHR_PREFIX and REFERENCE_GENOME

numberOfChromosomes=${CHROMOSOME_INDICES[@]}
outputFilenamePrefix=${mpileupDirectory}/${SNVFILE_PREFIX}${PID}

if [[ "$GERMLINE_AVAILABLE" == 0 ]]; then
    FILTER_VALUES=""
    [[ ${FILTER_ExAC} == 'true' ]]          && FILTER_VALUES="${FILTER_VALUES} ${ExAC_COL} AF ${CRIT_ExAC_maxMAF}+"
    [[ ${FILTER_EVS} == 'true' ]]           && FILTER_VALUES="${FILTER_VALUES} ${EVS_COL} MAF ${CRIT_EVS_maxMAF}+"
    #[[ ${FILTER_1KGENOMES} == 'true' ]]     && FILTER_VALUES="${FILTER_VALUES} ${KGENOMES_COL} AF ${CRIT_1KGENOMES_maxMAF}+"
    #[[ ${FILTER_1KGENOMES} == 'true' ]]     && FILTER_VALUES="${FILTER_VALUES} ${KGENOMES_COL} ASN_AF ${CRIT_1KGENOMES_maxMAF}+"
    #[[ ${FILTER_1KGENOMES} == 'true' ]]     && FILTER_VALUES="${FILTER_VALUES} ${KGENOMES_COL} AMR_AF ${CRIT_1KGENOMES_maxMAF}+"
    #[[ ${FILTER_1KGENOMES} == 'true' ]]     && FILTER_VALUES="${FILTER_VALUES} ${KGENOMES_COL} AFR_AF ${CRIT_1KGENOMES_maxMAF}+"
    [[ ${FILTER_1KGENOMES} == 'true' ]]     && FILTER_VALUES="${FILTER_VALUES} ${KGENOMES_COL} EUR_AF ${CRIT_1KGENOMES_maxMAF}+"
    [[ ${FILTER_NON_CLINIC} == 'true' ]]    && FILTER_VALUES="${FILTER_VALUES} ${DBSNP_COL} CLN,COMMON nonexist,exist"
    [[ ${FILTER_LOCALCONTROL} == 'true' ]]  && FILTER_VALUES="${FILTER_VALUES} ${LOCALCONTROL_COL} AF ${CRIT_LOCALCONTROL_maxMAF}+"
    [[ ${FILTER_RECURRENCE} == 'true' ]]    && FILTER_VALUES="${FILTER_VALUES} ${RECURRENCE_COL} . ${CRIT_RECURRENCE}+"

    if [[ ${FILTER_VALUES} != "" ]]; then
        echo ${FILTER_VALUES} > "${outputFilenamePrefix}_postFilter_criteria.txt"
        outputFilenamePrefix="${outputFilenamePrefix}_postFiltered"
        FILTERED_VCF="${outputFilenamePrefix}.vcf"
        ${PYPY_BINARY} -u ${TOOL_VCF_FILTER_BY_CRIT} ${FILENAME_VCF} ${FILTERED_VCF}${FILTER_VALUES}
        ${BGZIP_BINARY} -f ${FILTERED_VCF} && ${TABIX_BINARY} -f -p vcf ${FILTERED_VCF}.gz
        FILENAME_VCF=${FILTERED_VCF}.gz
    fi
fi

# file paths
filenameSomaticSnvs=${outputFilenamePrefix}_somatic_snvs_conf_${MIN_CONFIDENCE_SCORE}_to_10.vcf
filenameSomaticSnvsIndbSNP=${outputFilenamePrefix}_somatic_in_dbSNP_conf_${MIN_CONFIDENCE_SCORE}_to_10.txt
filenameIntermutationDistance=${outputFilenamePrefix}_somatic_mutation_dist_conf_${MIN_CONFIDENCE_SCORE}_to_10.txt.tmp
filenamePCRerrorMatrix=${outputFilenamePrefix}_sequence_specific_error_Matrix_conf_${MIN_CONFIDENCE_SCORE}_to_10.txt
filenameSequencingErrorMatrix=${outputFilenamePrefix}_sequencing_specific_error_Matrix_conf_${MIN_CONFIDENCE_SCORE}_to_10.txt
filenameReferenceAlleleBaseQualities=${outputFilenamePrefix}_reference_allele_base_qualities.txt.gz
filenameAlternativeAlleleBaseQualities=${outputFilenamePrefix}_alternative_allele_base_qualities.txt.gz
filenameAlternativeAlleleReadPositions=${outputFilenamePrefix}_alternative_allele_read_positions.txt.gz

# plot paths
filenamePerChromFreq=${outputFilenamePrefix}_perChromFreq_conf_${MIN_CONFIDENCE_SCORE}_to_10.pdf
filenameSnvsWithContext=${outputFilenamePrefix}_snvs_with_context_conf_${MIN_CONFIDENCE_SCORE}_to_10.pdf
filenameIntermutationDistancePlot=${outputFilenamePrefix}_intermutation_distance_conf_${MIN_CONFIDENCE_SCORE}_to_10.pdf
filenameSequenceErrorPlot=${outputFilenamePrefix}_sequence_specific_error_plot_conf_${MIN_CONFIDENCE_SCORE}_to_10.pdf
filenameSequencingErrorPlot=${outputFilenamePrefix}_sequencing_specific_error_plot_conf_${MIN_CONFIDENCE_SCORE}_to_10.pdf
filenameSequenceErrorPlotPreFilter=${outputFilenamePrefix}_sequence_specific_error_plot_before_filter.pdf
filenameSequencingErrorPlotPreFilter=${outputFilenamePrefix}_sequencing_specific_error_plot_before_filter.pdf
filenameSequenceErrorPlotFilterOnce=${outputFilenamePrefix}_sequence_specific_error_plot_after_filter_once.pdf
filenameSequencingErrorPlotFilterOnce=${outputFilenamePrefix}_sequencing_specific_error_plot_after_filter_once.pdf
filenameBaseScoreBiasPlotPreFilter=${outputFilenamePrefix}_base_score_bias_before_filter.pdf
filenameBaseScoreBiasPlotOnce=${outputFilenamePrefix}_base_score_bias_after_filter_once.pdf
filenameBaseScoreBiasPlotFinal=${outputFilenamePrefix}_base_score_bias_plot_conf_${MIN_CONFIDENCE_SCORE}_to_10.pdf
filenameBaseScoreDistributions=${outputFilenamePrefix}_base_score_distribution.pdf

# maf plots
filenameMafValues=${outputFilenamePrefix}_MAF_conf_${MIN_CONFIDENCE_SCORE}_to_10.txt.tmp
filenameSeqContextTab=${outputFilenamePrefix}_snvs_with_context_conf_${MIN_CONFIDENCE_SCORE}_to_10.txt.tmp
filenameMAFconfPlot=${outputFilenamePrefix}_MAF_conf_${MIN_CONFIDENCE_SCORE}_to_10.pdf
filenameSnvDiagnosticsPlot=${outputFilenamePrefix}_allSNVdiagnosticsPlots.pdf

${PERL_BINARY} ${TOOL_SNV_EXTRACTOR} --infile=${FILENAME_VCF} --minconf=${MIN_CONFIDENCE_SCORE} --pid=${outputFilenamePrefix} --bgzip=${BGZIP_BINARY} --tabix=${TABIX_BINARY} ${SNV_FILTER_OPTIONS}
[[ "$?" != 0 ]] && echo "There was a non-zero exit code in the somatic file and dbSNP counting pipe" && exit 1

if [ ${RUN_PLOTS} == 1 ]
then
	# get clinical info on PID => decoded PID, gender, age, subgroup

    # Test if CLINICALANNO is set to more than ""
	[[ -n "$CLINICALANNO" ]] && anno=`grep -w $PID $CLINICALANNO | awk '{print $2"_"$3"_"$4"_"$5}'` || anno=""
	[[ -n "$CLINICALANNO" ]] && gender=`grep -w $PID $CLINICALANNO | cut -f3 -s` || gender=""
	exY=''
	if [ "$gender" = F ] || [ "$gender" = f ]
	then
		ignoreY=1
		exY="--excludedChromosomes=chrY"
	fi

	# 1. rainfall plots
	cat ${filenameSomaticSnvs} | ${PERL_BINARY} ${TOOL_IN_DB_SNP_COUNTER} - ${MIN_CONFIDENCE_SCORE} > ${filenameSomaticSnvsIndbSNP}

	[[ "$?" != 0 ]] && echo "There was a non-zero exit code in the somatic file and dbSNP counting pipe" && exit 2

	# mutation distance, rainfall plot,  mutation classes per chromosome
    ${PYTHON_BINARY} ${TOOL_MUTATION_DISTANCE} --inf=${filenameSomaticSnvs} --outf=${filenameIntermutationDistance} --alleleFreq=${ALLELE_FREQ} ${exY}
    [[ "$?" != 0 ]] && echo "There was a non-zero exit code in making the mutation distance file" && exit 3
	
	CHR_TO_RAINFALL=${numberOfChromosomes}

	${RSCRIPT_BINARY} --vanilla ${TOOL_INTERMUTATION_DISTANCE_COORD_COLOR} -i ${filenameSomaticSnvs} -s ${PID} -o ${filenameIntermutationDistancePlot} -a ${CHR_TO_RAINFALL// /,} -p "${CHR_PREFIX}" -u "${CHR_SUFFIX}" -l ${CHROMOSOME_LENGTH_FILE}
	
	[[ "$?" != 0 ]] && echo "There was a non-zero exit code in making the rainfall plot file" && exit 4
	
	${RSCRIPT_BINARY} ${TOOL_SNVS_PER_CHROM_PLOT} -i ${filenameIntermutationDistance} -s ${PID}_${anno} -o ${filenamePerChromFreq}

	[[ "$?" != 0 ]] && echo "There was a non-zero exit code in making the SNVs per chromosome plot file" && exit 5

	# 2. MAF plots
	# get mutant allele frequency ("MAF") for extracted somatic high confidence SNVs
	${PERL_BINARY} ${TOOL_MAKE_MAF_INPUT} ${filenameSomaticSnvs} "$MINCOV" "$MIN_CONFIDENCE_SCORE" > $filenameMafValues

	[[ "$?" != 0 ]] && echo "There was a non-zero exit code in making the MAF input file" && exit 6

	# count the obtained SNVs to output their number in the plot: < 50 will not be reliable!
	snvnum=`grep -v "^#" ${filenameSomaticSnvs} | wc -l`
	snvindbSNP=` awk '{FS="\t"}{if(NR==2)print $5}'	${filenameSomaticSnvsIndbSNP}`

	# make MAF plot - from Natalie
	if [ "$snvnum" != "0" ]; then
		${RSCRIPT_BINARY} ${TOOL_MAF_PLOTS} ${filenameMafValues} ${snvnum} ${filenameMAFconfPlot} ${PID} ${snvindbSNP}

		[[ "$?" != 0 ]] && echo "There was a non-zero exit code in MAF plotting" && exit 7
	else
		filenameMAFconfPlot="" # no output produced, don't include in "convert" later
	fi

	# generate sequence contexts with preceding and following base for SNVs - from Matthias. Needs the original SNV file
	${PERL_BINARY} ${TOOL_SNV_CONTEXT_FREQUENCIES} ${filenameSomaticSnvs} $MIN_CONFIDENCE_SCORE > $filenameSeqContextTab
	[[ "$?" != 0 ]] && echo "There was a non-zero exit code in making the sequence context file" && exit 8

	# make sequence context plot - from Matthias Schlesner
	${RSCRIPT_BINARY} ${TOOL_SNV_SEQ_CONTEXT} ${filenameSeqContextTab} somatic ${filenameSnvsWithContext}

	[[ "$?" != 0 ]] && echo "There was a non-zero exit code in making the snv context plot" && exit 9

    rm ${filenameMafValues} ${filenameSeqContextTab} ${filenameIntermutationDistance}

    # make error plots from Matthias Bieg
    ${PYTHON_BINARY} ${TOOL_CREATE_ERROR_PLOTS} --vcfFile=${filenameSomaticSnvs} --referenceFile=NA --outputFile=${filenameSequencingErrorPlot} --errorType=sequencing_specific --errorFile=${filenameSequencingErrorMatrix} --plot_title='Final sequencing strand bias from vcf filter script'
    
    [[ "$?" != 0 ]] && echo "There was a non-zero exit code in making the sequencing bias files" && exit 10
    
    ${PYTHON_BINARY} ${TOOL_CREATE_ERROR_PLOTS} --vcfFile=${filenameSomaticSnvs} --referenceFile=NA --outputFile=${filenameSequenceErrorPlot} --errorType=sequence_specific --errorFile=${filenamePCRerrorMatrix} --plot_title='Final PCR strand bias from vcf filter script'

	[[ "$?" != 0 ]] && echo "There was a non-zero exit code in making the PCR bias files" && exit 11

    # make base score bias and base score distribution plots - by Gregor Warsow
    if [[ -f ${filenameReferenceAlleleBaseQualities} ]] && [[ -f ${filenameAlternativeAlleleBaseQualities} ]]; then
        basequal=`echo ${MPILEUP_OPTS} | perl -ne '($qual) = $_ =~ /\-Q\s*(\d+)/;print $qual'`
        # really use base score threshold of 13 if variable is not set pr empty? We should expect to get a proper value here...
        basequal=${basequal:-13}

#        ${RSCRIPT_BINARY} ${TOOL_PLOT_BASE_SCORE_BIAS} -v ${filenameSomaticSnvs} -r ${filenameReferenceAlleleBaseQualities} -a ${filenameAlternativeAlleleBaseQualities} -o ${filenameBaseScoreBiasPlot} -d "Final Base Quality Bias Plot for PID ${PID}"
        ${RSCRIPT_BINARY} ${TOOL_PLOT_BASE_SCORE_BIAS} -v ${filenameSomaticSnvs} -r ${filenameReferenceAlleleBaseQualities} -a ${filenameAlternativeAlleleBaseQualities} -t ${basequal} -p ${PLOT_TYPE} -o ${filenameBaseScoreBiasPlotFinal} -d "Final Base Quality Bias Plot for PID ${PID}"

        ${RSCRIPT_BINARY} ${TOOL_PLOT_BASE_SCORE_DISTRIBUTION} -v ${filenameSomaticSnvs} -r ${filenameReferenceAlleleBaseQualities} -a ${filenameAlternativeAlleleBaseQualities} -o ${filenameBaseScoreDistributions} -d "for somatic SNVs for PID ${PID}" -t ${basequal}
    else
        filenameBaseScoreDistributions=''
        filenameBaseScoreBiasPlotPreFilter=''
        filenameBaseScoreBiasPlotOnce=''
        filenameBaseScoreBiasPlotFinal=''
    fi

	# make a pdf containing all plots
	biasplots=""
	[[ -f ${filenameSequencingErrorPlot} ]] && biasplots="${filenameSequencingErrorPlot} ${biasplots}"
	[[ -f ${filenameSequencingErrorPlotFilterOnce} ]] && biasplots="${filenameSequencingErrorPlotFilterOnce} ${biasplots}"
	[[ -f ${filenameSequencingErrorPlotPreFilter} ]] && biasplots="${filenameSequencingErrorPlotPreFilter} ${biasplots}"
	[[ -f ${filenameSequenceErrorPlot} ]] && biasplots="${filenameSequenceErrorPlot} ${biasplots}"
	[[ -f ${filenameSequenceErrorPlotFilterOnce} ]] && biasplots="${filenameSequenceErrorPlotFilterOnce} ${biasplots}"
	[[ -f ${filenameSequenceErrorPlotPreFilter} ]] && biasplots="${filenameSequenceErrorPlotPreFilter} ${biasplots}"
	[[ -f ${filenameBaseScoreBiasPlotPreFilter} ]] && biasplots="${biasplots} ${filenameBaseScoreBiasPlotPreFilter}"
	[[ -f ${filenameBaseScoreBiasPlotOnce} ]] && biasplots="${biasplots} ${filenameBaseScoreBiasPlotOnce}"
	[[ -f ${filenameBaseScoreBiasPlotFinal} ]] && biasplots="${biasplots} ${filenameBaseScoreBiasPlotFinal}"
	
	${GHOSTSCRIPT_BINARY} -dBATCH -dNOPAUSE -dAutoRotatePages=false -q -sDEVICE=pdfwrite -sOutputFile=${filenameSnvDiagnosticsPlot} ${filenameIntermutationDistancePlot} ${filenameBaseScoreDistributions} ${filenameMAFconfPlot} ${filenamePerChromFreq} ${filenameSnvsWithContext} ${biasplots}
fi

if [ ${RUN_PUREST} == 1 ]
then
	# 3. purityEST - from Florian. Needs the original SNV file because it also considers germline (DP5 field)
	# has everything hardcoded (in which fields to look and confidence 8)
	confCol=`${PERL_BINARY} ${TOOL_FIND_CONF_COLUMN} ${FILENAME_VCF}`
	${PYTHON_BINARY} ${TOOL_PURITY_RELOADED} ${FILENAME_VCF} ${confCol} > ${outputFilenamePrefix}_purityEST.txt
	[[ "$?" != 0 ]] && echo "There was a non-zero exit code in purity estimation" && exit 7
fi

touch ${FILENAME_CHECKPOINT}
