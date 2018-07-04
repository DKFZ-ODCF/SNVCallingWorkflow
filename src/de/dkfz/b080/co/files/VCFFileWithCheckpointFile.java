package de.dkfz.b080.co.files;

import de.dkfz.roddy.knowledge.files.BaseFile;
import de.dkfz.roddy.knowledge.files.FileGroup;

import java.util.List;

/**
 * A vcf file with a checkpoint file. This is to simulate a tuple for generic method calls.
 */
public class VCFFileWithCheckpointFile extends FileGroup<BaseFile> {

    private SNVAnnotationFile vcfFile;
    private TextFile checkpointFile;

    public VCFFileWithCheckpointFile(List<BaseFile> files) {
        super(files);
        for (int i = 0; i < files.size(); i++) {
            BaseFile file = files.get(i);
            if (file instanceof SNVAnnotationFile)
                vcfFile = (SNVAnnotationFile) file;
            else if (file instanceof TextFile)
                checkpointFile = (TextFile) file;
        }
    }

    public SNVAnnotationFile getVCFFile() {
        return vcfFile;
    }

    public TextFile getCheckpointFile() {
        return checkpointFile;
    }
}
