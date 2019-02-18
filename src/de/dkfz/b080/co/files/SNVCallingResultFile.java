/*
 * Copyright (c) 2018 German Cancer Research Center (DKFZ).
 *
 * Distributed under the MIT License (license terms are at https://github.com/DKFZ-ODCF/IndelCallingWorkflow).
 *
 */
package de.dkfz.b080.co.files;

import de.dkfz.roddy.knowledge.files.BaseFile;

/**
 *
 * @author michael
 */
public abstract class SNVCallingResultFile extends BaseFile {

    public SNVCallingResultFile(ConstructionHelperForBaseFiles helper) {
        super(helper);
    }
}
