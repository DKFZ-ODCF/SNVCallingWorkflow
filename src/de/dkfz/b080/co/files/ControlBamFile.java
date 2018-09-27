/*
 * Copyright (c) 2018 German Cancer Research Center (DKFZ).
 *
 * Distributed under the MIT License (license terms are at https://github.com/DKFZ-ODCF/IndelCallingWorkflow).
 *
 */
package de.dkfz.b080.co.files;

import de.dkfz.roddy.knowledge.files.BaseFile;

/**
 * Created by heinold on 15.01.16.
 */
public class ControlBamFile extends BasicBamFile {

    public ControlBamFile(ConstructionHelperForBaseFiles helper) {
        super(helper);
    }

    /** "Copy" constructor */
    public ControlBamFile(BaseFile parent) {
        super(parent);
    }
}
