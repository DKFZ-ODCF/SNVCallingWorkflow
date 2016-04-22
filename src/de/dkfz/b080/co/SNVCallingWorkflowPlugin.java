package de.dkfz.b080.co;

import de.dkfz.roddy.plugins.BasePlugin;

/**
 */
public class SNVCallingWorkflowPlugin extends BasePlugin {

    public static final String CURRENT_VERSION_STRING = "1.1.2";
    public static final String CURRENT_VERSION_BUILD_DATE = "Fri Apr 22 15:36:08 CEST 2016";

    @Override
    public String getVersionInfo() {
        return "Roddy plugin: " + this.getClass().getName() + ", V " + CURRENT_VERSION_STRING + " built at " + CURRENT_VERSION_BUILD_DATE;
    }
}
