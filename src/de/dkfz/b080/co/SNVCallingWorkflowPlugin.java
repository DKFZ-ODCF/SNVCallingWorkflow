package de.dkfz.b080.co;

import de.dkfz.roddy.plugins.BasePlugin;

/**
 */
public class SNVCallingWorkflowPlugin extends BasePlugin {

    public static final String CURRENT_VERSION_STRING = "1.3.2";
    public static final String CURRENT_VERSION_BUILD_DATE = "Fri Nov 10 09:34:33 CET 2017";

    @Override
    public String getVersionInfo() {
        return "Roddy plugin: " + this.getClass().getName() + ", V " + CURRENT_VERSION_STRING + " built at " + CURRENT_VERSION_BUILD_DATE;
    }
}
