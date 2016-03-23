package de.dkfz.b080.co.snvpipeline;

import de.dkfz.b080.co.common.WorkflowUsingMergedBams;
import de.dkfz.b080.co.files.*;
import de.dkfz.roddy.core.ExecutionContext;

/**
 */
public class SNVCallingWorkflow extends WorkflowUsingMergedBams {

    @Override
    public boolean execute(ExecutionContext context, BasicBamFile _bamControlMerged, BasicBamFile _bamTumorMerged) {

        BamFile bamControlMerged = new BamFile(_bamControlMerged);
        BamFile bamTumorMerged = new BamFile(_bamTumorMerged);

        boolean runMetaCallingStep = context.getConfiguration().getConfigurationValues().getBoolean("runSNVMetaCallingStep", false);
        boolean runDeepAnnotation = context.getConfiguration().getConfigurationValues().getBoolean("runDeepAnnotation", true);
        boolean runFilter = context.getConfiguration().getConfigurationValues().getBoolean("runFilter", true);

        SNVAnnotationFile rawVCFFile = null;
        if(!runMetaCallingStep) {
            VCFFileGroupForSNVs vcfFilesForSNVs = Methods.callSNVs(bamControlMerged, bamTumorMerged);
            rawVCFFile = vcfFilesForSNVs.join();
        } else {
            rawVCFFile = Methods.callSNVsMeta(bamControlMerged, bamTumorMerged);
        }

        VCFFileWithCheckpointFile annotationFile = rawVCFFile.annotate();
        if (runDeepAnnotation)
            annotationFile = annotationFile.getVCFFile().deepAnnotate();
        if(runFilter)
            annotationFile = annotationFile.getVCFFile().filter(rawVCFFile, bamTumorMerged);

        return true;
    }
}
