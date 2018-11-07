/*
* To change this template, choose Tools | Templates
* and open the template in the editor.
*/
package de.dkfz.b080.co.files;

import de.dkfz.roddy.knowledge.files.FileGroup;

import java.util.Arrays;
import java.util.List;

/**
* @author michael
*/
public class SNVCallingFileGroup extends FileGroup<SNVCallingResultFile> {

    private final VCFFileForSNVs vcfFileForSNVs;
//    private final List<MPileupErrorFile> errorFiles;

    public SNVCallingFileGroup(VCFFileForSNVs vcfFileForSNVs) {
        super(Arrays.<SNVCallingResultFile>asList(vcfFileForSNVs));
        this.vcfFileForSNVs = vcfFileForSNVs;
//        this.vcfFileForIndels = vcfFileForIndels;
    }

    public SNVCallingFileGroup(List<SNVCallingResultFile> files) {
        super(files);
        //TODO Can there be other file input combinations?
        this.vcfFileForSNVs = (VCFFileForSNVs)(files.get(0) instanceof VCFFileForSNVs ? files.get(0) : files.get(1));
//        this.vcfFileForIndels = (VCFFileForIndels)(files.get(1) instanceof VCFFileForIndels ? files.get(1) : files.get(0));
    }

    public VCFFileForSNVs getVCFFileForSNVs() {
        return vcfFileForSNVs;
    }
//
//    public VCFFileForIndels getVCFFileForIndels() {
//        return vcfFileForIndels;
//    }
//
    @Override
    public void runDefaultOperations() {
//        deepAnnotateSNVFiles();
//        deepAnnotateIndelFiles();
    }

}
