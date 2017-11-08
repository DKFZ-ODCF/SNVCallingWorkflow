package de.dkfz.b080.co;

import de.dkfz.roddy.plugins.BasePlugin;

/**
 */
public class SNVCallingWorkflowPlugin extends BasePlugin {

    public static final String CURRENT_VERSION_STRING = "1.2.166";
    public static final String CURRENT_VERSION_BUILD_DATE = "Fri Aug 18 19:36:35 CEST 2017";

    @Override
    public String getVersionInfo() {
        return "Roddy plugin: " + this.getClass().getName() + ", V " + CURRENT_VERSION_STRING + " built at " + CURRENT_VERSION_BUILD_DATE;
    }
}
