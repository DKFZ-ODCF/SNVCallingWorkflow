package de.dkfz.b080.co.snvpipeline;

/**
 * Created by warsow on 4/18/18.
 */
@groovy.transform.CompileStatic
public class SNVCallingConstants {

    /**
     * Tool entries
     */
    public static final String TOOL_SNV_CALLING = "snvCalling";
    public static final String TOOL_SNV_CALLING_NOCONTROL = "snvCallingNoControl";
    public static final String TOOL_JOIN_SNV_VCF_FILES = "snvJoinVcfFiles";
    public static final String TOOL_SNV_ANNOTATION = "snvAnnotation";
    public static final String TOOL_SNV_DEEP_ANNOTATION = "snvDeepAnnotation";
    public static final String TOOL_SNV_FILTER = "snvFilter";
    public static final String TOOL_SNV_FILTER_RERUN = "snvFilterRerun";


    public static final String PARM_CHR_INDEX = "PARM_CHR_INDEX";
    public static final String CHR_NAME = "CHR_NAME";
    public static final String CHR_NR = "CHR_NR";
//    public static final String GENETIC_MAP_FILE = "GENETIC_MAP_FILE";
//    public static final String KNOWN_HAPLOTYPES_FILE = "KNOWN_HAPLOTYPES_FILE";
//    public static final String KNOWN_HAPLOTYPES_LEGEND_FILE = "KNOWN_HAPLOTYPES_LEGEND_FILE";

    private SNVCallingConstants() {
    }
}
