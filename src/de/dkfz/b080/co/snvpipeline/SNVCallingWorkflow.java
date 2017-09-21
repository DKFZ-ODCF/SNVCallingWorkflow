package de.dkfz.b080.co.snvpipeline;

import de.dkfz.b080.co.common.WorkflowUsingMergedBams;
import de.dkfz.b080.co.files.*;
import de.dkfz.roddy.core.ExecutionContext;
import de.dkfz.roddy.knowledge.files.Tuple2;

/**
 */
public class SNVCallingWorkflow extends WorkflowUsingMergedBams {

    @Override
    public boolean execute(ExecutionContext context, BasicBamFile _bamControlMerged, BasicBamFile _bamTumorMerged) {

        ControlBamFile bamControlMerged = new ControlBamFile(_bamControlMerged);
        TumorBamFile bamTumorMerged = new TumorBamFile(_bamTumorMerged);

        boolean runMetaCallingStep =  context.getConfiguration().getConfigurationValues().getBoolean("runSNVMetaCallingStep", false);
        boolean runDeepAnnotation = context.getConfiguration().getConfigurationValues().getBoolean("runDeepAnnotation", true);
        boolean runFilter = context.getConfiguration().getConfigurationValues().getBoolean("runFilter", true);
        final boolean runSecondFilterStep = getflag(context,"runSecondFilterStep", false);

        SNVAnnotationFile rawVCFFile = null;
        if(!runMetaCallingStep) {
            VCFFileGroupForSNVs vcfFilesForSNVs = Methods.callSNVs(bamControlMerged, bamTumorMerged);
            rawVCFFile = vcfFilesForSNVs.join();
        } else {
            rawVCFFile = (SNVAnnotationFile) call("snvCallingMetaScript", bamTumorMerged, bamControlMerged);
        }

        Tuple2<SNVAnnotationFile,TextFile> firstFilterRunResult = null;
        VCFFileWithCheckpointFile annotationFile = rawVCFFile.annotate();
        if (runDeepAnnotation)
            annotationFile = annotationFile.getVCFFile().deepAnnotate();
        if(runFilter)
            firstFilterRunResult = annotationFile.getVCFFile().filter(rawVCFFile, bamTumorMerged);
        if (firstFilterRunResult != null && runSecondFilterStep)
            firstFilterRunResult.value0.filterRerun(rawVCFFile, bamTumorMerged, firstFilterRunResult.value1);

        return true;
    }
}
