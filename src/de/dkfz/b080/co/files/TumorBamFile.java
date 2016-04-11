package de.dkfz.b080.co.files;

import de.dkfz.roddy.knowledge.files.BaseFile;

/**
 * Created by heinold on 15.01.16.
 */
public class TumorBamFile extends BasicBamFile {

    public TumorBamFile(ConstructionHelperForBaseFiles helper) {
        super(helper);
    }

    /** "Copy" constructor */
    public TumorBamFile(BaseFile parent) {
        super(parent);
    }
}
