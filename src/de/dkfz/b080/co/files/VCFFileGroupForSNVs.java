/*
* To change this template, choose Tools | Templates
* and open the template in the editor.
*/
package de.dkfz.b080.co.files;

import de.dkfz.roddy.knowledge.files.FileGroup;
import de.dkfz.roddy.knowledge.methods.GenericMethod;

import java.util.List;

//import pipelines.snv.SNVTools;

/**
 * @author michael
 */
public class VCFFileGroupForSNVs extends FileGroup<VCFFileForSNVs> {

    public VCFFileGroupForSNVs(List<VCFFileForSNVs> files) {
        super(files);
    }

    public SNVAnnotationFile join() {
        BamFile bf = (BamFile) filesInGroup.get(0).getParentFiles().get(0); //Should be merged tumor bam file
        SNVAnnotationFile file = GenericMethod.callGenericTool(COConstants.TOOL_JOIN_SNV_VCF_FILES, bf, this);
        return file;
//    }
//
//    @Override
//    public SNVAnnotationFile annotate() {
//        SNVAnnotationFile file = GenericMethod.callGenericTool(COConstants.TOOL_SNV_ANNOTATION, bf, this);
//        return file;
//    }
//
//    @Override
//    public SNVAnnotationFile deepAnnotate() {
//    }

    }
}
