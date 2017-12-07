package de.dkfz.b080.co.snvpipeline;

import de.dkfz.b080.co.common.WorkflowUsingMergedBams;
import de.dkfz.b080.co.files.BamFile;
import de.dkfz.b080.co.files.SNVAnnotationFile;
import de.dkfz.b080.co.files.VCFFileGroupForSNVs;
import de.dkfz.b080.co.files.VCFFileWithCheckpointFile;
import de.dkfz.roddy.config.ConfigurationError;
import de.dkfz.roddy.core.ExecutionContext;

/**
 */
public class SNVCallingWorkflow extends WorkflowUsingMergedBams {

    @Override
    public boolean execute(ExecutionContext context, BamFile bamControlMerged, BamFile bamTumorMerged) throws ConfigurationError {

        boolean runMetaCallingStep = context.getConfiguration().getConfigurationValues().getBoolean("runSNVMetaCallingStep", false);
        boolean runDeepAnnotation = context.getConfiguration().getConfigurationValues().getBoolean("runDeepAnnotation", true);
        boolean runFilter = context.getConfiguration().getConfigurationValues().getBoolean("runFilter", true);

        SNVAnnotationFile rawVCFFile = null;
        if(!runMetaCallingStep) {
            VCFFileGroupForSNVs vcfFilesForSNVs = bamTumorMerged.callSNVs(bamControlMerged);
            rawVCFFile = vcfFilesForSNVs.join();
        } else {
            rawVCFFile = bamTumorMerged.callSNVsMeta(bamControlMerged);
        }
        VCFFileWithCheckpointFile annotationFile = rawVCFFile.annotate();
        if (runDeepAnnotation)
            annotationFile = annotationFile.getVCFFile().deepAnnotate();
        if(runFilter)
            annotationFile = annotationFile.getVCFFile().filter(rawVCFFile, bamTumorMerged);

        return true;
    }
}
