package de.dkfz.b080.co.snvpipeline;

import de.dkfz.b080.co.common.COConfig;
import de.dkfz.b080.co.common.WorkflowUsingMergedBams;
import de.dkfz.b080.co.files.*;
import de.dkfz.roddy.core.ExecutionContext;
import de.dkfz.roddy.knowledge.files.Tuple2;
import de.dkfz.roddy.tools.LoggerWrapper;

/**
 */
public class SNVCallingWorkflow extends WorkflowUsingMergedBams {

private static LoggerWrapper logger = LoggerWrapper.getLogger(SNVCallingWorkflow.class.getName());
    @Override
    public boolean execute(ExecutionContext context, BasicBamFile _bamControlMerged, BasicBamFile _bamTumorMerged) {
        logger.postAlwaysInfo("noControlFLAG is " + isNoControlWorkflow());
        BasicBamFile bamTumorMerged;
        if (isControlWorkflow()) {
            BasicBamFile bamControlMerged = new BasicBamFile(loadInitialBamFilesForDataset(context)[0]);
            context.getConfigurationValues().add(new ConfigurationValue("controlSample", ((COFileStageSettings) bamControlMerged.getFileStage()).getSample().getName()));
            bamTumorMerged = new BasicBamFile(loadInitialBamFilesForDataset(context)[1]);
        } else {
            bamTumorMerged = new BasicBamFile(loadInitialBamFilesForDataset(context)[0]);
        }
        context.getConfigurationValues().add(new ConfigurationValue("tumorSample", ((COFileStageSettings) bamTumorMerged.getFileStage()).getSample().getName()));


        boolean runMetaCallingStep = context.getConfiguration().getConfigurationValues().getBoolean("runSNVMetaCallingStep", false);
        boolean runDeepAnnotation = context.getConfiguration().getConfigurationValues().getBoolean("runDeepAnnotation", true);
        boolean runFilter = context.getConfiguration().getConfigurationValues().getBoolean("runFilter", true);
        final boolean runSecondFilterStep = context.getConfiguration().getConfigurationValues().getBoolean("runSecondFilterStep", false);

        SNVAnnotationFile rawVCFFile = null;
        if (runMetaCallingStep && !noControlFLAG) {
            rawVCFFile = (SNVAnnotationFile) call("snvCallingMetaScript", bamTumorMerged, bamControlMerged);
        } else {
            final VCFFileGroupForSNVs vcfFilesForSNVs;
            if (noControlFLAG) {
                vcfFilesForSNVs = Methods.callSNVsNoControl(bamTumorMerged);
            } else {
                vcfFilesForSNVs = Methods.callSNVs(bamControlMerged, bamTumorMerged);
            }
            rawVCFFile = vcfFilesForSNVs.join();
        }

        Tuple2<SNVAnnotationFile, TextFile> firstFilterRunResult = null;
        VCFFileWithCheckpointFile annotationFile = rawVCFFile.annotate();
        if (runDeepAnnotation)
            annotationFile = annotationFile.getVCFFile().deepAnnotate();
        if (runFilter)
            firstFilterRunResult = annotationFile.getVCFFile().filter(rawVCFFile, bamTumorMerged);
        if (firstFilterRunResult != null && runSecondFilterStep)
            firstFilterRunResult.value0.filterRerun(rawVCFFile, bamTumorMerged, firstFilterRunResult.value1);

        return true;
    }
}
