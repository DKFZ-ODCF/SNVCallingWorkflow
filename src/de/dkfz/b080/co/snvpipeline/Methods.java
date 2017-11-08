package de.dkfz.b080.co.snvpipeline;

import de.dkfz.b080.co.common.ParallelizationHelper;
import de.dkfz.b080.co.files.*;
import de.dkfz.roddy.config.Configuration;
import de.dkfz.roddy.execution.jobs.ScriptCallingMethod;
import de.dkfz.roddy.knowledge.files.IndexedFileObjects;
import de.dkfz.roddy.knowledge.files.Tuple2;
import de.dkfz.roddy.knowledge.methods.GenericMethod;

import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;

/**
 * Created by heinold on 13.01.16.
 */
public class Methods {

    final String CHR_NAME = "CHR_NAME";
    final String CHR_NR = "CHR_NR";

    static LinkedHashMap<String, String> getGlobalJobSpecificParameters(Configuration config) {
        LinkedHashMap<String, String> results = new LinkedHashMap<String, String>();
//        results.put(CHR_NAME, config.getConfigurationValues().getString(CHR_NAME));
//        results.put(CHR_NR, config.getConfigurationValues().getString(CHR_NR));
        return results;
    }

    @ScriptCallingMethod
    public static VCFFileGroupForSNVs callSNVsNoControl(TumorBamFile tumorBam) {
        IndexedFileObjects<Tuple2<VCFFileForSNVs, TextFile>> indexedFileObjects = ParallelizationHelper.runParallel(COConstants.CVALUE_CHROMOSOME_INDICES, "snvCallingNoControl", tumorBam, null, "PARM_CHR_INDEX", getGlobalJobSpecificParameters(tumorBam.getExecutionContext().getConfiguration()));

        return getVcfFileGroupForSNVs(indexedFileObjects);
    }

    @ScriptCallingMethod
    public static VCFFileGroupForSNVs callSNVs(ControlBamFile controlBam, TumorBamFile tumorBam) {
        IndexedFileObjects<Tuple2<VCFFileForSNVs, TextFile>> indexedFileObjects = ParallelizationHelper.runParallel(COConstants.CVALUE_CHROMOSOME_INDICES, COConstants.TOOL_SNV_CALLING, tumorBam, controlBam, "PARM_CHR_INDEX", getGlobalJobSpecificParameters(controlBam.getExecutionContext().getConfiguration()));

        return getVcfFileGroupForSNVs(indexedFileObjects);
    }

    private static VCFFileGroupForSNVs getVcfFileGroupForSNVs(IndexedFileObjects<Tuple2<VCFFileForSNVs, TextFile>> indexedFileObjects) {
        for (String chromosome : indexedFileObjects.getIndices()) {
            indexedFileObjects.getIndexedFileObjects().get(chromosome);
        }

        List<VCFFileForSNVs> vcfFileForSNVsList = new LinkedList<>();
        for (String chromosomeIndex : indexedFileObjects.getIndices()) {
            Tuple2<VCFFileForSNVs, TextFile> tuple = indexedFileObjects.getIndexedFileObjects().get(chromosomeIndex);
            VCFFileForSNVs vcfFile = tuple.value0;
            TextFile checkpointFile = tuple.value1;
            vcfFile.setAsTemporaryFile();
            vcfFileForSNVsList.add(vcfFile);
        }
        VCFFileGroupForSNVs vcfFilesForSNVs = new VCFFileGroupForSNVs(vcfFileForSNVsList);
        return vcfFilesForSNVs;
    }

}
