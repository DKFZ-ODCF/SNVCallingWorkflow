#!/bin/bash

source ${CONFIG_FILE}

PLOT_TYPE=${BASE_SCORE_BIAS_PLOT_TYPE:-Differences}

if [[ ${FILENAME_CHECKPOINT_FIRST_FILTER_RUN:-0} == 0 ]]; then
    RERUN_FILTER_STEP=0
    MEDIAN_FILTER_THRESHOLD=-1
    RERUN_SUFFIX=""
else
    RERUN_FILTER_STEP=1
fi


set -o pipefail

[[ -f ${FILENAME_CHECKPOINT} ]] && rm ${FILENAME_CHECKPOINT}
mpileupDirectory=`dirname ${FILENAME_VCF}`

RUN_PLOTS=${RUN_PLOTS-1}
RUN_PUREST=${RUN_PUREST-1}
CONVERT_BINARY=${CONVERT_BINARY-convert}

source ${TOOL_ANALYZE_BAM_HEADER}
getRefGenomeAndChrPrefixFromHeader ${TUMOR_BAMFILE_FULLPATH_BP} # Sets CHR_PREFIX and REFERENCE_GENOME
ALIGNMENT_FOLDER=`dirname ${TUMOR_BAMFILE_FULLPATH_BP}`

# bugfix: ensure to interpret CHROMOSOME_INDICES as array - otherwise TOOL_INTERMUTATION_DISTANCE_COORD_COLOR will fail...
declare -a CHROMOSOME_INDICES="${CHROMOSOME_INDICES}"
numberOfChromosomes=${CHROMOSOME_INDICES[@]}
outputFilenamePrefix=${mpileupDirectory}/${SNVFILE_PREFIX}${PID}
outputFilenamePrefix_original=${outputFilenamePrefix}
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
filenameSomaticSnvs=${outputFilenamePrefix}_somatic_snvs_conf_${MIN_CONFIDENCE_SCORE}_to_10${RERUN_SUFFIX}.vcf
filenameSomaticSnvsIndbSNP=${outputFilenamePrefix}_somatic_in_dbSNP_conf_${MIN_CONFIDENCE_SCORE}_to_10${RERUN_SUFFIX}.txt
filenameIntermutationDistance=${outputFilenamePrefix}_somatic_mutation_dist_conf_${MIN_CONFIDENCE_SCORE}_to_10.txt.tmp
filenamePCRerrorMatrix=${outputFilenamePrefix}_sequence_specific_error_Matrix_conf_${MIN_CONFIDENCE_SCORE}_to_10${RERUN_SUFFIX}.txt
filenameSequencingErrorMatrix=${outputFilenamePrefix}_sequencing_specific_error_Matrix_conf_${MIN_CONFIDENCE_SCORE}_to_10${RERUN_SUFFIX}.txt
filenameReferenceAlleleBaseQualities=${outputFilenamePrefix_original}_reference_allele_base_qualities.txt.gz
filenameReferenceAlleleReadPositions=${outputFilenamePrefix_original}_reference_allele_read_positions.txt.gz
filenameAlternativeAlleleBaseQualities=${outputFilenamePrefix_original}_alternative_allele_base_qualities.txt.gz
filenameAlternativeAlleleReadPositions=${outputFilenamePrefix_original}_alternative_allele_read_positions.txt.gz
filenameTHAArtifactDetected=${outputFilenamePrefix}_is_THA_affected.txt
filenameQCvalues=${outputFilenamePrefix}_QC_values${RERUN_SUFFIX}.json

# plot paths
filenamePerChromFreq=${outputFilenamePrefix}_perChromFreq_conf_${MIN_CONFIDENCE_SCORE}_to_10${RERUN_SUFFIX}.pdf
filenameSnvsWithContext=${outputFilenamePrefix}_snvs_with_context_conf_${MIN_CONFIDENCE_SCORE}_to_10${RERUN_SUFFIX}.pdf
filenameIntermutationDistancePlot=${outputFilenamePrefix}_intermutation_distance_conf_${MIN_CONFIDENCE_SCORE}_to_10${RERUN_SUFFIX}.pdf
filenameSequenceErrorPlot=${outputFilenamePrefix}_sequence_specific_error_plot_conf_${MIN_CONFIDENCE_SCORE}_to_10${RERUN_SUFFIX}.pdf
filenameSequencingErrorPlot=${outputFilenamePrefix}_sequencing_specific_error_plot_conf_${MIN_CONFIDENCE_SCORE}_to_10${RERUN_SUFFIX}.pdf
filenameSequenceErrorPlotPreFilter=${outputFilenamePrefix_original}_sequence_specific_error_plot_before_filter.pdf
filenameSequencingErrorPlotPreFilter=${outputFilenamePrefix_original}_sequencing_specific_error_plot_before_filter.pdf
filenameSequenceErrorPlotFilterOnce=${outputFilenamePrefix_original}_sequence_specific_error_plot_after_filter_once.pdf
filenameSequencingErrorPlotFilterOnce=${outputFilenamePrefix_original}_sequencing_specific_error_plot_after_filter_once.pdf
filenameBaseScoreBiasPlotPreFilter=${outputFilenamePrefix_original}_base_score_bias_before_filter.pdf
filenameBaseScoreBiasPlotOnce=${outputFilenamePrefix_original}_base_score_bias_after_filter_once.pdf
filenameBaseScoreBiasPlotFinal=${outputFilenamePrefix}_base_score_bias_plot_conf_${MIN_CONFIDENCE_SCORE}_to_10${RERUN_SUFFIX}.pdf
filenameBaseScoreDistributions=${outputFilenamePrefix}_base_score_distribution${RERUN_SUFFIX}.pdf
BaseScoreDistributionsPlots_PREFIX=${outputFilenamePrefix}_tripletSpecific_base_score_distribution${RERUN_SUFFIX}
BaseScoreDistributionsPlots_COMBINED=${BaseScoreDistributionsPlots_PREFIX}.BQD_combined.pdf
BaseScoreDistributionsPlots_CoV=${BaseScoreDistributionsPlots_PREFIX}.BQD_CoV.pdf
BaseScoreDistributionsPlots_INDIVIDUAL=${BaseScoreDistributionsPlots_PREFIX}.BQD_individual.pdf
BaseScoreDistributionsPlots_INDIVIDUAL_ChromColored=${BaseScoreDistributionsPlots_PREFIX}.BQD_individual_CHROMcolored.pdf
BaseScoreDistributionsPlots_INDIVIDUAL_VAFColored=${BaseScoreDistributionsPlots_PREFIX}.BQD_individual_VAFcolored.pdf
BaseScoreDistributionsPlots_INDIVIDUAL_ReadPositionColored=${BaseScoreDistributionsPlots_PREFIX}.BQD_individual_ReadPosColored_Q60.pdf

# maf plots
filenameMafValues=${outputFilenamePrefix}_MAF_conf_${MIN_CONFIDENCE_SCORE}_to_10.txt.tmp
filenameSeqContextTab=${outputFilenamePrefix}_snvs_with_context_conf_${MIN_CONFIDENCE_SCORE}_to_10.txt.tmp
filenameMAFconfPlot=${outputFilenamePrefix}_MAF_conf_${MIN_CONFIDENCE_SCORE}_to_10${RERUN_SUFFIX}.pdf
filenameSnvDiagnosticsPlot=${outputFilenamePrefix}_allSNVdiagnosticsPlots${RERUN_SUFFIX}.pdf



if [[ ${RERUN_FILTER_STEP} == 1 ]]; then

    # prepare call of base score dist plot script in order to get median-filtered vcf file (skip plotting itself)
    plotBackgroundBaseScoreDistribution='0'
    forceRerun='1'
    combineRevComp='1'
    channelIndividualGraphs='1'
    skipPlots='1'

    filenameSomaticSnvs_original=`basename ${outputFilenamePrefix}_somatic_snvs_conf_${MIN_CONFIDENCE_SCORE}_to_10.vcf`
    SEQUENCE_CONTEXT_COLUMN_INDEX=`cat ${mpileupDirectory}/${filenameSomaticSnvs_original} | grep -v '^##' | grep '^#' | perl -ne 'use List::Util qw(first); chomp; my @colnames = split(/\t/, $_); my $columnIndex = first { $colnames[$_] eq "SEQUENCE_CONTEXT"} 0..$#colnames; $columnIndex += 1; print "$columnIndex\n";'`
    SNV_FILE_WITH_MAF=${mpileupDirectory}/${filenameSomaticSnvs_original}.withMAF.vcf
    SNV_FILE_WITH_MAF_filtered=${mpileupDirectory}/${filenameSomaticSnvs_original}.withMAF_filteredAltMedian${MEDIAN_FILTER_THRESHOLD}.vcf
    if [[ ! -f ${SNV_FILE_WITH_MAF} ]]; then
        cat ${mpileupDirectory}/${filenameSomaticSnvs_original} | perl -ne 'chomp; my $line=$_; if (/DP4=(\d+),(\d+),(\d+),(\d+);/) {my $fR=$1; my $rR=$2; my $fA=$3; my $rA=$4; my $MAF=($fA+$rA)/($fR+$rR+$fA+$rA); print "$line\t$MAF\n";} else { if (/^#CHROM/) { print "$line\tMAF\n";} else {print "$line\n";} };' >${SNV_FILE_WITH_MAF}
    fi
    MAF_COLUMN_INDEX=`cat ${SNV_FILE_WITH_MAF} | grep -v '^##' | grep '^#' | perl -ne 'use List::Util qw(first); chomp; my @colnames = split(/\t/, $_); my $columnIndex = first { $colnames[$_] eq "MAF"} 0..$#colnames; $columnIndex += 1; print "$columnIndex\n";'`
    # this script will create a file named ${SNV_FILE_WITH_MAF_filtered}
    ${RSCRIPT_BINARY} ${TOOL_PLOT_TRIPLET_SPECIFIC_BASE_SCORE_DISTRIBUTION} -v ${SNV_FILE_WITH_MAF} -m ${mpileupDirectory} -a ${ALIGNMENT_FOLDER} \
        -p ${PID} -b ${plotBackgroundBaseScoreDistribution} -o ${BaseScoreDistributionsPlots_PREFIX} -R ${forceRerun} -c ${combineRevComp} \
        -f ${MEDIAN_FILTER_THRESHOLD} -s ${SEQUENCE_CONTEXT_COLUMN_INDEX} --MAFColumnIndex ${MAF_COLUMN_INDEX} -i ${channelIndividualGraphs} \
        -t 'Base score distribution of PID '${PID}'\nafter Median'${MEDIAN_FILTER_THRESHOLD}' filtering' --skipPlots ${skipPlots} \
        --refBaseQual ${filenameReferenceAlleleBaseQualities} --altBaseQual ${filenameAlternativeAlleleBaseQualities} \
        --altReadPos ${filenameAlternativeAlleleReadPositions} --refReadPos ${filenameReferenceAlleleReadPositions}
    [[ "$?" != 0 ]] && "There was a non-zero exit code in the generation of the tripletBased BQ distribution plot." && exit 21
    mv ${SNV_FILE_WITH_MAF_filtered} ${filenameSomaticSnvs}

    cp ${filenameSomaticSnvs} ${filenameSomaticSnvs}.forSNVExtractor
    ${PERL_BINARY} ${TOOL_SNV_EXTRACTOR} --infile=${filenameSomaticSnvs}.forSNVExtractor --minconf=${MIN_CONFIDENCE_SCORE} --pid=${outputFilenamePrefix} --suffix=${RERUN_SUFFIX} --bgzip=${BGZIP_BINARY} --tabix=${TABIX_BINARY} ${SNV_FILTER_OPTIONS}
    [[ "$?" != 0 ]] && echo "There was a non-zero exit code in the somatic file and dbSNP counting pipe" && exit 25
    rm ${filenameSomaticSnvs}.forSNVExtractor
    rm ${SNV_FILE_WITH_MAF}

    filenameSomaticFunctionalSnvs_original=`basename ${outputFilenamePrefix}_somatic_functional_snvs_conf_${MIN_CONFIDENCE_SCORE}_to_10.vcf`
    filenameSomaticFunctionalSnvs_MedianFiltered=`basename ${outputFilenamePrefix}_somatic_functional_snvs_conf_${MIN_CONFIDENCE_SCORE}_to_10${RERUN_SUFFIX}.vcf`
    filenameSomaticFunctionalSnvs_RemovedByMedianFilter=`basename ${outputFilenamePrefix}_somatic_functional_snvs_conf_${MIN_CONFIDENCE_SCORE}_to_10_removedByMedian${MEDIAN_FILTER_THRESHOLD}Filter.vcf`
    grep '^#' ${filenameSomaticFunctionalSnvs_original} >${filenameSomaticFunctionalSnvs_RemovedByMedianFilter}; ${BEDTOOLS_BINARY} subtract -a ${filenameSomaticFunctionalSnvs_original} -b ${filenameSomaticFunctionalSnvs_MedianFiltered} >>${filenameSomaticFunctionalSnvs_RemovedByMedianFilter}
else
    ${PERL_BINARY} ${TOOL_SNV_EXTRACTOR} --infile=${FILENAME_VCF} --minconf=${MIN_CONFIDENCE_SCORE} --pid=${outputFilenamePrefix} --bgzip=${BGZIP_BINARY} --tabix=${TABIX_BINARY} ${SNV_FILTER_OPTIONS}
    [[ "$?" != 0 ]] && echo "There was a non-zero exit code in the somatic file and dbSNP counting pipe" && exit 26
fi


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
	
	${RSCRIPT_BINARY} ${TOOL_SNVS_PER_CHROM_PLOT} -i ${filenameIntermutationDistance} -l ${CHROMOSOME_LENGTH_FILE} -s ${PID}_${anno} -o ${filenamePerChromFreq}

	[[ "$?" != 0 ]] && echo "There was a non-zero exit code in making the SNVs per chromosome plot file" && exit 5

	# 2. MAF plots
	# get mutant allele frequency ("MAF") for extracted somatic high confidence SNVs
	${PERL_BINARY} ${TOOL_MAKE_MAF_INPUT} ${filenameSomaticSnvs} "$MINCOV" "$MIN_CONFIDENCE_SCORE" > $filenameMafValues

	[[ "$?" != 0 ]] && echo "There was a non-zero exit code in making the MAF input file" && exit 6

	# count the obtained SNVs to output their number in the plot: < 50 will not be reliable!
	snvnum=`grep -v "^#" ${filenameSomaticSnvs} | wc -l`
	snvindbSNP=` awk '{FS="\t"}{if(NR==2)print $5}'	${filenameSomaticSnvsIndbSNP}`
	# QC value $SNV_IN_DBSNP_RATIO will be written to $filenameQCvalues
    SNV_IN_DBSNP_RATIO=`echo -e "$snvindbSNP\t$snvnum" | perl -F -ne 'print $F[0]/$F[1];'`

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
        [[ "$?" != 0 ]] && echo "There was a non-zero exit code in the base score bias plot script." && exit 31
        if [[ ! -f ${filenameBaseScoreBiasPlotFinal} ]]; then
            filenameBaseScoreBiasPlotFinal=''
        fi

        ${RSCRIPT_BINARY} ${TOOL_PLOT_BASE_SCORE_DISTRIBUTION} -v ${filenameSomaticSnvs} -r ${filenameReferenceAlleleBaseQualities} -a ${filenameAlternativeAlleleBaseQualities} -o ${filenameBaseScoreDistributions} -d "for somatic SNVs for PID ${PID}" -t ${basequal}
        [[ "$?" != 0 ]] && echo "There was a non-zero exit code in the base score distribution plots script." && exit 32
        if [[ ! -f ${filenameBaseScoreDistributions} ]]; then
            filenameBaseScoreDistributions=''
        fi

        plotBackgroundBaseScoreDistribution='0'
        if [[ ${RERUN_FILTER_STEP} == 1 ]]; then
            # do not force rerun as data file has been created in upper part (get median-filtered vcf file)
            forceRerun='0'
        else
            forceRerun='1'
        fi
        combineRevComp='1'
        channelIndividualGraphs='1'
        skipPlots='0'
        SNV_FILE_WITH_MAF=${filenameSomaticSnvs}.withMAF.vcf
        cat ${filenameSomaticSnvs} | perl -ne 'chomp; my $line=$_; if (/DP4=(\d+),(\d+),(\d+),(\d+);/) {my $fR=$1; my $rR=$2; my $fA=$3; my $rA=$4; my $MAF=($fA+$rA)/($fR+$rR+$fA+$rA); print "$line\t$MAF\n";} else { if (/^#CHROM/) { print "$line\tMAF\n";} else {print "$line\n";} };' >${SNV_FILE_WITH_MAF}

        SEQUENCE_CONTEXT_COLUMN_INDEX=`cat ${SNV_FILE_WITH_MAF} | grep -v '^##' | grep '^#' | perl -ne 'use List::Util qw(first); chomp; my @colnames = split(/\t/, $_); my $columnIndex = first { $colnames[$_] eq "SEQUENCE_CONTEXT"} 0..$#colnames; $columnIndex += 1; print "$columnIndex\n";'`
        MAF_COLUMN_INDEX=`cat ${SNV_FILE_WITH_MAF} | grep -v '^##' | grep '^#' | perl -ne 'use List::Util qw(first); chomp; my @colnames = split(/\t/, $_); my $columnIndex = first { $colnames[$_] eq "MAF"} 0..$#colnames; $columnIndex += 1; print "$columnIndex\n";'`
        ${RSCRIPT_BINARY} ${TOOL_PLOT_TRIPLET_SPECIFIC_BASE_SCORE_DISTRIBUTION} -v ${SNV_FILE_WITH_MAF} -m ${mpileupDirectory} -a ${ALIGNMENT_FOLDER} -p ${PID} \
        -b ${plotBackgroundBaseScoreDistribution} -o ${BaseScoreDistributionsPlots_PREFIX} -R ${forceRerun} -c ${combineRevComp} -f ${MEDIAN_FILTER_THRESHOLD} \
        -s ${SEQUENCE_CONTEXT_COLUMN_INDEX} --MAFColumnIndex ${MAF_COLUMN_INDEX} -i ${channelIndividualGraphs} -t 'Base score distribution of PID '${PID} \
        --skipPlots ${skipPlots} --refBaseQual ${filenameReferenceAlleleBaseQualities} --altBaseQual ${filenameAlternativeAlleleBaseQualities} \
        --altReadPos ${filenameAlternativeAlleleReadPositions} --refReadPos ${filenameReferenceAlleleReadPositions}
        if [[ ! ${runSecondFilterStep} ]]; then
            rm ${SNV_FILE_WITH_MAF}
        fi


        [[ ! -f ${BaseScoreDistributionsPlots_COMBINED} ]] && BaseScoreDistributionsPlots_COMBINED=''
        [[ ! -f ${BaseScoreDistributionsPlots_CoV} ]] && BaseScoreDistributionsPlots_CoV=''
        [[ ! -f ${BaseScoreDistributionsPlots_INDIVIDUAL} ]] && BaseScoreDistributionsPlots_INDIVIDUAL=''
        [[ ! -f ${BaseScoreDistributionsPlots_INDIVIDUAL_ChromColored} ]] && BaseScoreDistributionsPlots_INDIVIDUAL_ChromColored=''
        [[ ! -f ${BaseScoreDistributionsPlots_INDIVIDUAL_VAFColored} ]] && BaseScoreDistributionsPlots_INDIVIDUAL_VAFColored=''
        [[ ! -f ${BaseScoreDistributionsPlots_INDIVIDUAL_ReadPositionColored} ]] && BaseScoreDistributionsPlots_INDIVIDUAL_ReadPositionColored=''

    else
        filenameBaseScoreDistributions=''
        filenameBaseScoreBiasPlotPreFilter=''
        filenameBaseScoreBiasPlotOnce=''
        filenameBaseScoreBiasPlotFinal=''
        BaseScoreDistributionsPlots_COMBINED=''
        BaseScoreDistributionsPlots_CoV=''
        BaseScoreDistributionsPlots_INDIVIDUAL=''
        BaseScoreDistributionsPlots_INDIVIDUAL_ChromColored=''
        BaseScoreDistributionsPlots_INDIVIDUAL_VAFColored=''
        BaseScoreDistributionsPlots_INDIVIDUAL_ReadPositionColored=''
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

	BaseScoreDistributionsPlots_INDIVIDUAL_ChromColored='' # turn off ChromColored BSD
	${GHOSTSCRIPT_BINARY} -dBATCH -dNOPAUSE -dAutoRotatePages=false -q -sDEVICE=pdfwrite -sOutputFile=${filenameSnvDiagnosticsPlot} ${filenameIntermutationDistancePlot} ${filenameMAFconfPlot} ${filenamePerChromFreq} ${filenameSnvsWithContext} ${filenameBaseScoreDistributions} ${BaseScoreDistributionsPlots_COMBINED} ${BaseScoreDistributionsPlots_INDIVIDUAL} ${BaseScoreDistributionsPlots_INDIVIDUAL_ChromColored} ${BaseScoreDistributionsPlots_INDIVIDUAL_VAFColored} ${BaseScoreDistributionsPlots_INDIVIDUAL_ReadPositionColored} ${BaseScoreDistributionsPlots_CoV} ${biasplots}
fi

    # infer baseQuality bias (PV4)-related THA score (QC value)
    THA_SCORE=`${RSCRIPT_BINARY} ${TOOL_THA_DETECTOR} -i ${filenameSomaticSnvs}`
    [[ "$?" != 0 ]] && echo "There was a non-zero exit code in THA score determination script." && exit 24
    [[ -f ${filenameTHAArtifactDetected} ]] && rm ${filenameTHAArtifactDetected}
    [[ $(echo "${THA_SCORE} > ${THA_SCORE_THRESHOLD}" | bc -l) ]] && echo -e "THA score\t${THA_SCORE}\n" >${filenameTHAArtifactDetected}


if [ ${RUN_PUREST} == 1 ]
then
	# 3. purityEST - from Florian. Needs the original SNV file because it also considers germline (DP5 field)
	# has everything hardcoded (in which fields to look and confidence 8)
	confCol=`${PERL_BINARY} ${TOOL_FIND_CONF_COLUMN} ${FILENAME_VCF}`
	${PYTHON_BINARY} ${TOOL_PURITY_RELOADED} ${FILENAME_VCF} ${confCol} > ${outputFilenamePrefix}_purityEST${RERUN_SUFFIX}.txt
	[[ "$?" != 0 ]] && echo "There was a non-zero exit code in purity estimation" && exit 7
fi

# determine fraction of SNVs called as "synonymous SNV" among all exonic SNVs (QC value)
EXONIC_CLASSIFICATION_COLUMN_INDEX=`cat ${filenameSomaticSnvs} | grep -v '^##' | grep '^#' | perl -ne 'use List::Util qw(first); chomp; my @colnames = split(/\t/, $_); my $columnIndex = first { $colnames[$_] eq "EXONIC_CLASSIFICATION"} 0..$#colnames; $columnIndex += 1; print "$columnIndex\n";'`
export EXONIC_CLASSIFICATION_COLUMN_INDEX=$((${EXONIC_CLASSIFICATION_COLUMN_INDEX}-1))
ANNOVAR_FUNCTION_COLUMN_INDEX=`cat ${filenameSomaticSnvs} | grep -v '^##' | grep '^#' | perl -ne 'use List::Util qw(first); chomp; my @colnames = split(/\t/, $_); my $columnIndex = first { $colnames[$_] eq "ANNOVAR_FUNCTION"} 0..$#colnames; $columnIndex += 1; print "$columnIndex\n";'`
export ANNOVAR_FUNCTION_COLUMN_INDEX=$((${ANNOVAR_FUNCTION_COLUMN_INDEX}-1))
SYNONYMOUS_RATIO=`grep -v '^#' ${filenameSomaticSnvs} | perl -F'\t' -ae 'BEGIN { my $total=0; my $synonymous=0; } if ($F[$ENV{"ANNOVAR_FUNCTION_COLUMN_INDEX"}] eq "exonic" ) {$total++; if ($F[$ENV{"EXONIC_CLASSIFICATION_COLUMN_INDEX"}] eq "synonymous SNV") {$synonymous++;}} END { print $synonymous/$total; }'`

echo -e "{" >$filenameQCvalues
echo -e "\t\"snvnum\": ${snvnum:-NA}," >>$filenameQCvalues
echo -e "\t\"snvindbSNP\": ${snvindbSNP:-NA}," >>$filenameQCvalues
echo -e "\t\"SNV_IN_DBSNP_RATIO\": ${SNV_IN_DBSNP_RATIO:-NA}," >>$filenameQCvalues
echo -e "\t\"SYNONYMOUS_RATIO\": ${SYNONYMOUS_RATIO:-NA}," >>$filenameQCvalues
echo -e "\t\"THA_SCORE\": ${THA_SCORE:-NA}" >>$filenameQCvalues
#echo -e "\t\"GOX_SCORE\": ${GOX_SCORE}" >>$filenameQCvalues
echo -e "}" >>$filenameQCvalues

touch ${FILENAME_CHECKPOINT}
