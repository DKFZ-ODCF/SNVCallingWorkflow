/*
 * Copyright (c) 2018 German Cancer Research Center (DKFZ).
 *
 * Distributed under the MIT License (license terms are at https://github.com/DKFZ-ODCF/IndelCallingWorkflow).
 *
 */
package de.dkfz.b080.co.files;

import de.dkfz.b080.co.snvpipeline.SNVCallingConstants;
import de.dkfz.roddy.knowledge.files.BaseFile;
import de.dkfz.roddy.knowledge.files.Tuple2;
import de.dkfz.roddy.knowledge.methods.GenericMethod;

/**
 * @author michael
 */
public class SNVAnnotationFile extends BaseFile {

    private BaseFile parentFile;

    public SNVAnnotationFile(ConstructionHelperForBaseFiles helper) {
        super(helper);
        if (helper instanceof ConstructionHelperForGenericCreation){
            this.parentFile = (BaseFile) ((ConstructionHelperForGenericCreation) helper).parentObject;}
    }

    public VCFFileWithCheckpointFile annotate() {
        VCFFileWithCheckpointFile file = GenericMethod.callGenericTool(SNVCallingConstants.TOOL_SNV_ANNOTATION, this, parentFile);
        return file;
    }

    public VCFFileWithCheckpointFile deepAnnotate() {
        VCFFileWithCheckpointFile file = GenericMethod.callGenericTool(SNVCallingConstants.TOOL_SNV_DEEP_ANNOTATION, this, "PIPENAME=SNV_DEEPANNOTATION");
        return file;
    }

    public Tuple2<SNVAnnotationFile,TextFile> filter(SNVAnnotationFile rawVCFFile, TumorBamFile tumorBamFile) {
        Tuple2<SNVAnnotationFile,TextFile> tuple = GenericMethod.callGenericTool(SNVCallingConstants.TOOL_SNV_FILTER, this, rawVCFFile, tumorBamFile);
        return tuple;
    }

    public Tuple2<SNVAnnotationFile,TextFile> filterRerun(SNVAnnotationFile rawVCFFile, TumorBamFile tumorBamFile, TextFile firstFilterRunCheckpointFile) {
        Tuple2<SNVAnnotationFile,TextFile> tuple = GenericMethod.callGenericTool(SNVCallingConstants.TOOL_SNV_FILTER_RERUN, this, rawVCFFile, tumorBamFile, firstFilterRunCheckpointFile);
        return tuple;
    }
}
