#!/bin/bash

set -o pipefail

[[ -f ${FILENAME_CHECKPOINT} ]] && rm ${FILENAME_CHECKPOINT}
mpileupDirectory=`dirname ${FILENAME_VCF}`

RUN_PLOTS=${RUN_PLOTS-1}
RUN_PUREST=${RUN_PUREST-1}

source ${TOOL_ANALYZE_BAM_HEADER}
getRefGenomeAndChrPrefixFromHeader ${TUMOR_BAMFILE_FULLPATH_BP} # Sets CHR_PREFIX and REFERENCE_GENOME

# bugfix: ensure to interpret CHROMOSOME_INDICES as array - otherwise TOOL_INTERMUTATION_DISTANCE_COORD_COLOR will fail...
declare -a CHROMOSOME_INDICES="${CHROMOSOME_INDICES}"
numberOfChromosomes=${CHROMOSOME_INDICES[@]}
outputFilenamePrefix=${mpileupDirectory}/${SNVFILE_PREFIX}${PID}

# file paths
filenameSomaticSnvs=${outputFilenamePrefix}_somatic_snvs_conf_${MIN_CONFIDENCE_SCORE}_to_10.vcf
filenameSomaticSnvsIndbSNP=${outputFilenamePrefix}_somatic_in_dbSNP_conf_${MIN_CONFIDENCE_SCORE}_to_10.txt
filenameIntermutationDistance=${outputFilenamePrefix}_somatic_mutation_dist_conf_${MIN_CONFIDENCE_SCORE}_to_10.txt.tmp
filenamePCRerrorMatrix=${outputFilenamePrefix}_sequence_specific_error_Matrix_conf_${MIN_CONFIDENCE_SCORE}_to_10.txt
filenameSequencingErrorMatrix=${outputFilenamePrefix}_sequencing_specific_error_Matrix_conf_${MIN_CONFIDENCE_SCORE}_to_10.txt

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

# maf plots
filenameMafValues=${outputFilenamePrefix}_MAF_conf_${MIN_CONFIDENCE_SCORE}_to_10.txt.tmp
filenameSeqContextTab=${outputFilenamePrefix}_snvs_with_context_conf_${MIN_CONFIDENCE_SCORE}_to_10.txt.tmp
filenameMAFconfPlot=${outputFilenamePrefix}_MAF_conf_${MIN_CONFIDENCE_SCORE}_to_10.pdf
filenameSnvDiagnosticsPlot=${outputFilenamePrefix}_allSNVdiagnosticsPlots.pdf

${PERL_BINARY} ${TOOL_SNV_EXTRACTOR} --infile=${FILENAME_VCF} --minconf=${MIN_CONFIDENCE_SCORE} --pid=${outputFilenamePrefix} --bgzip=${BGZIP_BINARY} --tabix=${TABIX_BINARY} --whitelist="${WHITELIST:-NA}" ${SNV_FILTER_OPTIONS}
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
	snvnum=`grep -v "#" ${filenameSomaticSnvs} | wc -l`
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

	# make a pdf containing all plots
	biasplots=""
	[[ -f ${filenameSequencingErrorPlot} ]] && biasplots="${filenameSequencingErrorPlot} ${biasplots}"
	[[ -f ${filenameSequencingErrorPlotFilterOnce} ]] && biasplots="${filenameSequencingErrorPlotFilterOnce} ${biasplots}"
	[[ -f ${filenameSequencingErrorPlotPreFilter} ]] && biasplots="${filenameSequencingErrorPlotPreFilter} ${biasplots}"
	[[ -f ${filenameSequenceErrorPlot} ]] && biasplots="${filenameSequenceErrorPlot} ${biasplots}"
	[[ -f ${filenameSequenceErrorPlotFilterOnce} ]] && biasplots="${filenameSequenceErrorPlotFilterOnce} ${biasplots}"
	[[ -f ${filenameSequenceErrorPlotPreFilter} ]] && biasplots="${filenameSequenceErrorPlotPreFilter} ${biasplots}"

	# Ghostscript does not want paths containing '+' characters (such as in timestamps).
	${GHOSTSCRIPT_BINARY} -dBATCH -dNOPAUSE -dAutoRotatePages=false -q -sDEVICE=pdfwrite -sOutputFile=- ${filenameIntermutationDistancePlot} ${filenameMAFconfPlot} ${filenamePerChromFreq} ${filenameSnvsWithContext} ${biasplots} \
	    > "$filenameSnvDiagnosticsPlot"
	[[ "$?" != 0 ]] && echo "There was a non-zero exit code in making the SNV diagnostic plot" && exit 12
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
