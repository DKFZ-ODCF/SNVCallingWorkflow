library(getopt)
library(ggplot2)
library(Biostrings) # for reverseComplement
library(grid) # for unit,gpar
library(reshape2)
library("RColorBrewer")


CANNEL_INDIVIDUAL_GRAPHS=T
MAX_BASE_QUALITY=55
ALT.MEDIAN.THRESHOLD = -1 # -1 means NO FILTERING
parameter.density.bandwidth = 2.00
SEQUENCE_CONTEXT_COLUMN_INDEX=11
VAF_COLUMN_INDEX=-1

opt = getopt(matrix(c(
  'vcfInputFile', 'v', 1, "character",
  'mpileupFolder', 'm', 1, "character",
  'alignmentFolder', 'a', 1, "character",
  'PID', 'p', 1, "character",
  'background', 'b', 2, "integer",
  'outFile', 'o', 1, "character",
  'forceRerun', 'R', 2, "integer",
  'combineRevcomp', 'c', 2, "integer",
  'filterThreshold', 'f', 2, "integer",
  'sequenceContextColumnIndex', 's', 2, "integer",
  'MAFColumnIndex', 'maf', 1, "integer",
  'channelIndividualGraphs', 'i', 2, "integer",
  'mainTitle', 't', 2, "character"
),ncol=4,byrow=TRUE));

if (is.null(opt$vcfInputFile)){
  cat("Please specify the file that contains the SNVs for which the base score distribution plot shall be created.\n"); 
  q(status=1);      # quit, status unequal 0 means error
}
if (is.null(opt$mpileupFolder)){
  cat("Please specify the mpileup folder.\n"); 
  q(status=1);      # quit, status unequal 0 means error
}
if (is.null(opt$alignmentFolder)){
  cat("Please specify the alignment folder.\n"); 
  q(status=1);      # quit, status unequal 0 means error
}
if (is.null(opt$PID)){
  cat("Please specify the PID.\n"); 
  q(status=1);      # quit, status unequal 0 means error
}
if (is.null(opt$background)){
  USE_BACKGROUND_BQs = FALSE
} else {
  if (opt$background==1) {
    USE_BACKGROUND_BQs = TRUE
  } else {
    USE_BACKGROUND_BQs = FALSE
  }
}
if (is.null(opt$forceRerun)){     
  FORCE_RERUN = FALSE
} else {
  if (opt$forceRerun==1) {
    FORCE_RERUN = TRUE
  } else {
    FORCE_RERUN = FALSE
  }
}
if (is.null(opt$combineRevcomp)){     
  COMBINE_REVCOMP = TRUE
} else {
  if (opt$combineRevcomp==1) {
    COMBINE_REVCOMP = TRUE
  } else {
    COMBINE_REVCOMP = FALSE
  }
}
if (! is.null(opt$MAFColumnIndex)){
  tmp = as.integer(opt$MAFColumnIndex)
  if (!is.na(tmp)) {
    VAF_COLUMN_INDEX = tmp
  }
}

if (is.null(opt$outFile)){      # no vcf file specified
  cat("Please specify the output pdf file.\n"); 
  q(status=2);      # quit, status unequal 0 means error
}

if (! is.null(opt$filterThreshold)){
  tmp = as.integer(opt$filterThreshold)
  if (!is.na(tmp)) {
    ALT.MEDIAN.THRESHOLD = tmp
  }
}
if (! is.null(opt$sequenceContextColumnIndex)){
  tmp = as.integer(opt$sequenceContextColumnIndex)
  if (!is.na(tmp)) {
    SEQUENCE_CONTEXT_COLUMN_INDEX = tmp
  }
}
if (! is.null(opt$channelIndividualGraphs)){
  if (opt$channelIndividualGraphs == 1) {
    CANNEL_INDIVIDUAL_GRAPHS = T
  } else {
    CANNEL_INDIVIDUAL_GRAPHS = F
  }
}




if (COMBINE_REVCOMP) {
  transitions = data.frame( c(rep("C",3),rep("T",3)), c("A","G","T","A","C","G"), stringsAsFactors = F)  
  possible_mutations = c("CA", "CG", "CT", "TA", "TC", "TG")
  forbidden_mutations = sapply(possible_mutations, function(mutation) {as.character(complement(DNAString(mutation)))})
} else {
  transitions = data.frame( c(rep("A",3),rep("C",3),rep("G",3),rep("T",3)), c("C","G","T","A","G","T","A","C","T","A","C","G"), stringsAsFactors = F)
  possible_mutations = c("AC","AG","AT", "CA","CG","CT", "GA","GC","GT", "TA","TC","TG")
}
colnames(transitions) = c("FROM", "TO")


# RPP_FOLDER="/icgc/pcawg/analysis/train2full/projects/CMDI-UK/"
# opt$PID = "f8e61a02-654c-c226-e040-11ac0d481b60"
# opt$mpileupFolder = paste0(RPP_FOLDER,opt$PID,"/merged_may2016Calls_filtered_BQD_plots")
# opt$vcfInputFile = paste0(opt$mpileupFolder,"/",opt$PID,"_somatic.snv_mnv.vcf_filtered.OnlyPassed.withoutHeader.vcf")
# opt$alignmentFolder = paste0(RPP_FOLDER,opt$PID,"/alignment")
# opt$outFile = paste0(opt$mpileupFolder,"/snvs_",opt$PID,"_triplet-specific_BQ_distributions_unfiltered.manual.pdf")
# opt$background = F
# opt$combineRevcomp = 1
# opt$filterThreshold = -1
# opt$MAFColumnIndex = 9
# opt$sequenceContextColumnIndex = 10



# RPP_FOLDER="/icgc/dkfzlsdf/analysis/B080/warsow/Cavathrombus/results_per_pid/BW/"
# opt$PID = "BW"
# opt$mpileupFolder = paste0(RPP_FOLDER,"SNV_Calling.backup/mpileup2.3_VRenalisAdjacent_vs_NormalKidney/")
# opt$vcfInputFile = paste0(opt$mpileupFolder,"/snvs_",opt$PID,"_somatic_snvs_conf_8_to_10.vcf")
# opt$alignmentFolder = paste0(RPP_FOLDER,"/alignment/")
# opt$outFile = paste0(opt$mpileupFolder,"snvs_",opt$PID,"_triplet-specific_BQ_distributions.pdf")
# opt$background = F
# opt$combineRevcomp = 1

# RPP_FOLDER="/icgc/pcawg/analysis/train2full/projects/BOCA-UK/"
# opt$PID = "f85ae2b7-cebf-17a2-e040-11ac0c48033a"
# opt$mpileupFolder = paste0(RPP_FOLDER,opt$PID,"/merged_may2016Calls_filtered_BQD_plots")
# opt$vcfInputFile = paste0(opt$mpileupFolder,"/",opt$PID,"_somatic.snv_mnv.vcf_filtered.OnlyPassed.vcf")
# opt$vcfInputFile = paste0(opt$mpileupFolder,"/",opt$PID,"_somatic.snv_mnv.vcf_filtered.OnlyPassed.withoutHeader.vcf")
# opt$outFile = paste0(opt$mpileupFolder,"/snvs_",opt$PID,"_SangerArtifcatDetected.txt")
# opt$combineRevcomp = 1


# opt$PID = "4116268"
# opt$PID = 4108101
# opt$PID = 4104893
# opt$PID = 4178345
# opt$baseFolder = "/icgc/dkfzlsdf/analysis/mmml/whole_genome_pcawg/results_per_pid"
# opt$vcfInputFile = paste0("snvs_",opt$PID,"_somatic_snvs_conf_8_to_10.vcf")
# opt$mpileupFolder = "mpileup"
# opt$outFile = paste0("snvs_",opt$PID,"_triplet-specific_BQ_distributions.pdf")
# opt$background = F

# BLCA-US (PCAWG) test config
# opt$PID = "8c619cbc-9e91-4716-9711-5236e55d8f46"
# opt$baseFolder = "/icgc/pcawg/analysis/train2full/projects/BLCA-US"
# opt$vcfInputFile = paste0("snvs_",opt$PID,"_somatic_snvs_conf_8_to_10.vcf")
# opt$vcfInputFile = "snvs_8c619cbc-9e91-4716-9711-5236e55d8f46_somatic_snvs_conf_8_to_10.vcf"
# opt$mpileupFolder = "BQD_plots"
# opt$outFile = paste0("snvs_",opt$PID,"_triplet-specific_BQ_distributions.pdf")
# opt$background = 0
# opt$forceRerun = 1

# BRCA-EU (PCAWG) test config
# opt$PID = "fc8130e0-0bfa-bba4-e040-11ac0c48328d"
# opt$baseFolder = "/icgc/pcawg/analysis/train2full/projects/BRCA-EU"
# opt$vcfInputFile = paste0("snvs_",opt$PID,"_somatic_snvs_conf_8_to_10.vcf")
# opt$mpileupFolder = "BQD_plots"
# opt$outFile = paste0("snvs_",opt$PID,"_triplet-specific_BQ_distributions.pdf")
# opt$background = T

# CLLE-ES (PCAWG) test config
# opt$PID = "19cd4360-8392-4bc2-ae88-fdc1335d886b"
# opt$baseFolder = "/icgc/pcawg/analysis/train2full/projects/CLLE-ES"
# opt$vcfInputFile = paste0("snvs_",opt$PID,"_somatic_snvs_conf_8_to_10.vcf")
# opt$mpileupFolder = "BQD_plots"
# opt$outFile = paste0("snvs_",opt$PID,"_triplet-specific_BQ_distributions.pdf")

# # LAML-US (PCAWG) test config
# opt$PID = "1193a9c4-5aab-4cd7-a690-60c96bd1172d"
# opt$baseFolder = "/icgc/pcawg/analysis/train2full/projects/LAML-US"
# opt$vcfInputFile = paste0("snvs_",opt$PID,"_somatic_snvs_conf_8_to_10.vcf")
# opt$mpileupFolder = "BQD_plots"
# opt$outFile = paste0("snvs_",opt$PID,"_triplet-specific_BQ_distributions.pdf")

# opt$vcfInputFile = paste0("snvs_",opt$PID,"_somatic_snvs_conf_8_to_10.vcf")
# opt$mpileupFolder = "BQD_plots"
# opt$outFile = paste0("snvs_",opt$PID,"_triplet-specific_BQ_distributions.pdf")
# (PCAWG) test config
# opt$PID = "320ebe6f-48bd-40f1-9be5-1401a3bff328"
# opt$baseFolder = "/icgc/pcawg/analysis/train2full/projects/LAML-US"
# opt$vcfInputFile = paste0("snvs_",opt$PID,"_somatic_snvs_conf_8_to_10.vcf")
# opt$mpileupFolder = "BQD_plots"
# opt$outFile = paste0("snvs_",opt$PID,"_triplet-specific_BQ_distributions.pdf")



# numberOfCores = 12
# library(doParallel)
# registerDoParallel(numberOfCores)

PID=opt$PID
# BASEFOLDER=paste0(opt$baseFolder,"/")
MPILEUP_FOLDER=paste0(opt$mpileupFolder,"/")
ALIGNMENT_FOLDER=paste0(opt$alignmentFolder,"/")
vcfInputFile = paste0(opt$vcfInputFile)
RefAlleleBaseQualitiesFile=paste0(MPILEUP_FOLDER,"snvs_",PID,"_reference_allele_base_qualities.txt.gz", collapse = "")
AltAlleleBaseQualitiesFile=paste0(MPILEUP_FOLDER,"snvs_",PID,"_alternative_allele_base_qualities.txt.gz", collapse = "")
PDF_OUTPUT_FILE=paste0(opt$outFile)
# PID="f8e61a02-654c-c226-e040-11ac0d481b60"
# MPILEUP_FOLDER="/icgc/pcawg/analysis/train2full/projects/CMDI-UK/f8e61a02-654c-c226-e040-11ac0d481b60/merged_may2016Calls_filtered_BQD_plots/"
# ALIGNMENT_FOLDER="/icgc/pcawg/analysis/train2full/projects/CMDI-UK/f8e61a02-654c-c226-e040-11ac0d481b60/alignment/"
# vcfInputFile = "/icgc/pcawg/analysis/train2full/projects/CMDI-UK/f8e61a02-654c-c226-e040-11ac0d481b60/merged_may2016Calls_filtered_BQD_plots/f8e61a02-654c-c226-e040-11ac0d481b60_somatic.snv_mnv.vcf_filtered.OnlyPassed.vcf"
# RefAlleleBaseQualitiesFile=paste0(MPILEUP_FOLDER,"snvs_",PID,"_reference_allele_base_qualities.txt.gz", collapse = "")
# AltAlleleBaseQualitiesFile=paste0(MPILEUP_FOLDER,"snvs_",PID,"_alternative_allele_base_qualities.txt.gz", collapse = "")
# PDF_OUTPUT_FILE=paste0(MPILEUP_FOLDER,"snvs_",PID,"_triplet-specific_BQ_distributions_filteredAltMedian20.pdf")


if (ALT.MEDIAN.THRESHOLD > -1) {
  DATA_RESULTS_FILE=paste0(MPILEUP_FOLDER,"BQ_TripletSpecific_MedianFiltered",ALT.MEDIAN.THRESHOLD,".RData")
  VCF_OUTPUT_FILE.filtered = sub(".vcf$", paste0("_filteredAltMedian",ALT.MEDIAN.THRESHOLD,".vcf"), vcfInputFile)  
} else {
  DATA_RESULTS_FILE=paste0(MPILEUP_FOLDER,"BQ_TripletSpecific.RData")
}


if ( file.exists(DATA_RESULTS_FILE) & ! FORCE_RERUN ) {
  load(file = DATA_RESULTS_FILE)
  if ( USE_BACKGROUND_BQs & !exists("data.bq.background.dens") ) {
    data.bq.background = read.table(paste0(ALIGNMENT_FOLDER,"BaseQualitiesForChosenPositions.1000000.1.txt.gz", collapse = ""))
    colnames(data.bq.background) = "BQ"
    data.bq.background.dens = density(data.bq.background$BQ, bw = parameter.density.bandwidth)
    remove(data.bq.background)
  }
} else {
  # GO THROUGH DATA PROCESSING STEPS (either due to missing data file (e.g. first run) or being asked to do it although data file is available)
  if (USE_BACKGROUND_BQs) {
    data.bq.background = read.table(paste0(ALIGNMENT_FOLDER,"BaseQualitiesForChosenPositions.1000000.1.txt.gz", collapse = ""))
    colnames(data.bq.background) = "BQ"
    data.bq.background.dens = density(data.bq.background$BQ, bw = parameter.density.bandwidth)
    remove(data.bq.background)
  }
  
  
  # READ IN DATA FOR CHANNEL-SPECIFIC BQs
  if (VAF_COLUMN_INDEX > -1) {
    if (VAF_COLUMN_INDEX > SEQUENCE_CONTEXT_COLUMN_INDEX) {
      #------------------------------------------------------------------------------------------------------------------------------------------------------------CHR----POS----REF---ALT---------------------BaseBefore-----BaseAfter---------------------VAF-----
      readFileCommand = paste0("cat ",vcfInputFile, " | grep -v '^##' | cut -f 1,2,4,5,",SEQUENCE_CONTEXT_COLUMN_INDEX,",",VAF_COLUMN_INDEX," | perl -ne '$_ =~ /^(.*?)\t(.*?)\t(.*?)\t(.*?)\t[ACGTNacgtn]{9}([ACGTNacgtn]),([ACGTNacgtn])[ACGTNacgtn]{9}\t([\\d\\.]*)/; print \"$1\t$2\t$3\t$4\t$5\t$6\t$7\n\";' ")
    } else {
      #---------------------------------------------------------------------------------------------------------------------------------------CHR----POS----REF---ALT----------------------------VAF------------------------BaseBefore-----BaseAfter
      readFileCommand = paste0("cat ",vcfInputFile, " | grep -v '^##' | cut -f 1,2,4,5,",VAF_COLUMN_INDEX,",",SEQUENCE_CONTEXT_COLUMN_INDEX," | perl -ne '$_ =~ /^(.*?)\t(.*?)\t(.*?)\t(.*?)\t([\\d\\.]*)\t[ACGTNacgtn]{9}([ACGTNacgtn]),([ACGTNacgtn])[ACGTNacgtn]{9}/; print \"$1\t$2\t$3\t$4\t$6\t$7\t$5\n\";' ")
    }
  } else {
    #---------------------------------------------------------------------------------------------------------------------------------------CHR----POS----REF---ALT---------------------BaseBefore-----BaseAfter
    readFileCommand = paste0("cat ",vcfInputFile, " | grep -v '^##' | cut -f 1,2,4,5,",SEQUENCE_CONTEXT_COLUMN_INDEX," | perl -ne '$_ =~ /^(.*?)\t(.*?)\t(.*?)\t(.*?)\t[ACGTNacgtn]{9}([ACGTNacgtn]),([ACGTNacgtn])[ACGTNacgtn]{9}/; print \"$1\t$2\t$3\t$4\t$5\t$6\n\";' ")
  }
  data.bq.triplet = read.table(pipe(readFileCommand), comment.char = '', sep = "\t", header = T, stringsAsFactors = F)
  if (VAF_COLUMN_INDEX > -1) {
    colnames(data.bq.triplet) = c("CHROM", "POS","REF","ALT","BaseBefore","BaseAfter","VAF")
  } else {
    colnames(data.bq.triplet) = c("CHROM", "POS","REF","ALT","BaseBefore","BaseAfter")
  }
  data.bq.triplet$Triplet = with(data.bq.triplet, paste0(BaseBefore,REF,ALT,BaseAfter))
  # VCF = data.bq.triplet
  
  data.bq.ref = read.table(file = RefAlleleBaseQualitiesFile, stringsAsFactors = F)
  data.bq.alt = read.table(file = AltAlleleBaseQualitiesFile, stringsAsFactors = F)
  colnames(data.bq.ref) = c("CHROM","POS","BQ_string")
  colnames(data.bq.alt) = c("CHROM","POS","BQ_string")
  
  data.bq.triplet = merge(data.bq.triplet, data.bq.ref, by=c("CHROM","POS"), all.x = T)
  colnames(data.bq.triplet)[ncol(data.bq.triplet)] = "BQ_string.ref"
  data.bq.triplet = merge(data.bq.triplet, data.bq.alt, by=c("CHROM","POS"), all.x = T)
  colnames(data.bq.triplet)[ncol(data.bq.triplet)] = "BQ_string.alt"
  data.bq.triplet = data.bq.triplet[order(data.bq.triplet$Triplet),]
  
  
  data.bq.triplet.processed = lapply(seq(nrow(data.bq.triplet)), function(i) {
    line = data.bq.triplet[i,]
    baseScores.ref = as.integer(unlist(strsplit(unlist(strsplit(as.character(line[1,"BQ_string.ref"]), ";"))[1], ",")))
    baseScores.alt = as.integer(unlist(strsplit(unlist(strsplit(as.character(line[1,"BQ_string.alt"]), ";"))[1], ","))) 
    
    countsLine.ref = data.frame(matrix(0, ncol = MAX_BASE_QUALITY+1, nrow = 1))
    colnames(countsLine.ref) = c(paste0("BQ.ref.",sapply(seq(MAX_BASE_QUALITY+1)-1, function(i) {i})))
    countsLine.alt = data.frame(matrix(0, ncol = MAX_BASE_QUALITY+1, nrow = 1))
    colnames(countsLine.alt) = c(paste0("BQ.alt.",sapply(seq(MAX_BASE_QUALITY+1)-1, function(i) {i})))    
    counts.ref = table(baseScores.ref);  if (length(counts.ref)>0) {countsLine.ref[1,paste0("BQ.ref.",names(counts.ref))] = as.integer(unlist(counts.ref))}
    counts.alt = table(baseScores.alt);  if (length(counts.alt)>0) {countsLine.alt[1,paste0("BQ.alt.",names(counts.alt))] = as.integer(unlist(counts.alt))}
    
    line[1,"nBQ.ref"] = length(baseScores.ref)
    line[1,"nBQ.alt"] = length(baseScores.alt)
    
    baseScores.ref.mean = mean(baseScores.ref)
    baseScores.alt.mean = mean(baseScores.alt)
    line[1,"MeanBQ.ref"] = baseScores.ref.mean
    line[1,"MeanBQ.alt"] = baseScores.alt.mean    
    
    baseScores.common.mean = mean(c(baseScores.ref, baseScores.alt))
    line[1,"MeanBQ.common"] = baseScores.common.mean
    
    baseScores.common.median = median(c(baseScores.ref, baseScores.alt))
    line[1,"MeadianBQ.common"] = baseScores.common.median
    
    baseScores.ref.median = median(baseScores.ref)
    baseScores.alt.median = median(baseScores.alt)
    line[1,"MeadianBQ.ref"] = baseScores.ref.median
    line[1,"MeadianBQ.alt"] = baseScores.alt.median
    
    baseScores.ref.CoV = sd(baseScores.ref)/mean(baseScores.ref)*100
    baseScores.alt.CoV = sd(baseScores.alt)/mean(baseScores.alt)*100
    line[1,"CoV.ref"] = baseScores.ref.CoV
    line[1,"CoV.alt"] = baseScores.alt.CoV
    
    baseScores.ref.AuC.BQ30 = sum(counts.ref[as.integer(names(counts.ref))<=30])/sum(counts.ref)
    baseScores.alt.AuC.BQ30 = sum(counts.alt[as.integer(names(counts.alt))<=30])/sum(counts.alt)
    baseScores.ref.AuC.BQ25 = sum(counts.ref[as.integer(names(counts.ref))<=25])/sum(counts.ref)
    baseScores.alt.AuC.BQ25 = sum(counts.alt[as.integer(names(counts.alt))<=25])/sum(counts.alt)
    baseScores.ref.AuC.BQ20 = sum(counts.ref[as.integer(names(counts.ref))<=20])/sum(counts.ref)
    baseScores.alt.AuC.BQ20 = sum(counts.alt[as.integer(names(counts.alt))<=20])/sum(counts.alt)
    line[1,"AuC.ref.BQ20"] = baseScores.ref.AuC.BQ20
    line[1,"AuC.alt.BQ20"] = baseScores.alt.AuC.BQ20
    line[1,"AuC.ref.BQ25"] = baseScores.ref.AuC.BQ25
    line[1,"AuC.alt.BQ25"] = baseScores.alt.AuC.BQ25
    line[1,"AuC.ref.BQ30"] = baseScores.ref.AuC.BQ30
    line[1,"AuC.alt.BQ30"] = baseScores.alt.AuC.BQ30
    
    line = cbind(line, countsLine.ref, countsLine.alt)
    
    return (line)
  })
  
  
  data.bq.triplet = do.call(rbind, data.bq.triplet.processed)
  
  remove(data.bq.triplet.processed, data.bq.ref, data.bq.alt)
  
  data.bq.triplet$"BaseBefore" = factor(data.bq.triplet$"BaseBefore", levels=c("A","C","G","T"))
  data.bq.triplet$"REF" = factor(data.bq.triplet$"REF", levels=c("A","C","G","T"))
  data.bq.triplet$"ALT" = factor(data.bq.triplet$"ALT", levels=c("A","C","G","T"))
  data.bq.triplet$"BaseAfter" = factor(data.bq.triplet$"BaseAfter", levels=c("A","C","G","T"))
  if ( ! is.null(data.bq.triplet$VAF) ) {
    data.bq.triplet$"VAF" = as.numeric(as.character(data.bq.triplet$"VAF"))
  }
  data.bq.triplet$"Triplet" = as.factor(as.character(data.bq.triplet$"Triplet"))
  
  if (USE_BACKGROUND_BQs) {
    save(data.bq.triplet, data.bq.background.dens, file = DATA_RESULTS_FILE)
    # write.table(data.bq.triplet.filtered)
  } else {
    save(data.bq.triplet, file = DATA_RESULTS_FILE)
  }
}


data.bq.triplet$REVCOMP_TRIPLET = apply(data.bq.triplet, 1, function(line) {
  if (COMBINE_REVCOMP) {
    if (line["REF"] == "C" | line["REF"] == "T") {
      return(line["Triplet"])
    } else {
      triplet = line["Triplet"]
      triplet.revComp.tmp = as.character(reverseComplement(DNAString(triplet)))
      triplet.revComp = paste0(substr(triplet.revComp.tmp,1,1),substr(triplet.revComp.tmp,3,3),substr(triplet.revComp.tmp,2,2),substr(triplet.revComp.tmp,4,4))
      return(triplet.revComp)
    }    
  } else {
    return(line["Triplet"])
  }
})


# Filter SNVs by median ALT BQ
if (ALT.MEDIAN.THRESHOLD > -1) {
  data.bq.triplet.filtered = do.call(rbind, apply(transitions, 1, function(TRANSITION) {
    if ( "THIS SNV IS ARTIFACTED" == "THIS SNV IS ARTIFACTED" ) {
      # TODO: only filter channels that show artifact (criterion still has to be defined)
      SNVs = data.bq.triplet[with(data.bq.triplet, substr(REVCOMP_TRIPLET,2,2)==TRANSITION[1] & substr(REVCOMP_TRIPLET,3,3)==TRANSITION[2]),]
      indices = which(SNVs$MeadianBQ.alt >= ALT.MEDIAN.THRESHOLD)
      SNVs = SNVs[indices,]    
      return (SNVs)
    }
  }))
}
# data.bq.triplet.filtered = data.bq.triplet.filtered[order(data.bq.triplet.filtered$MeadianBQ.alt),]
# any(is.na(data.bq.triplet.filtered$nBQ.ref))



######################################################
# PLOTTING PART BEGINS HERE
######################################################

## Add an alpha value to a colour
add.alpha = function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}


redColor = rgb(1,0,0)
greenColor = rgb(0,1,0)

xlim = 45
offset.row=2
offset.col=2
POSITIONS.BEFORE=c("A"=4+offset.row,"C"=3+offset.row,"G"=2+offset.row,"T"=1+offset.row)
POSITIONS.FROMTO.COL=c("CA"=0,"CG"=1,"CT"=2,"TA"=0,"TC"=1,"TG"=2)
POSITIONS.FROMTO.ROW=c("CA"=0,"CG"=0,"CT"=0,"TA"=1,"TC"=1,"TG"=1)
POSITIONS.AFTER=c("A"=1+offset.col,"C"=2+offset.col,"G"=3+offset.col,"T"=4+offset.col)


plot_BQD_to_pdf = function(PDF_OUTPUT_FILE, data.bq.triplet, whatToPlot=c("BQD","BQ_CoV","BQD_sampleIndividual")) {
  pdf(PDF_OUTPUT_FILE, 8.27, 6.0)
  # pdf(PDF_OUTPUT_FILE, paper = "a4r")
  
    if (CANNEL_INDIVIDUAL_GRAPHS == T) {
      if (COMBINE_REVCOMP == F) {
        stop("Channel-individual plots are only supported for combined revcomp!")
      }
      
      grid.newpage(recording = F)
      if (whatToPlot == "BQD") {
        gl <- grid.layout(nrow=13, ncol=17, widths = unit(c(3,1,2.6,2,2,2,1,2,2,2,2,1,2,2,2,2,3), "null"), heights = unit(c(2,2,2,2,2,2.5,2.5,2,2,2,2.5,0.5), "null"))
      } else if (whatToPlot == "BQ_CoV" || whatToPlot == "BQD_sampleIndividual") {
        gl <- grid.layout(nrow=13, ncol=17, widths = unit(c(3,1,2,2,2,2,1,2,2,2,2,1,2,2,2,2,3), "null"), heights = unit(c(2,2,2,2,2,2,2.5,2,2,2,2,0.5), "null"))
      }
      
      pushViewport(viewport(layout = gl))
      
      grid.text(opt$mainTitle, vp = viewport(layout.pos.row = 1, layout.pos.col = 2+offset.col), hjust = -0.2, vjust = 1.7)
      
      grid.text(paste0(transitions[1,1],"->",transitions[1,2]), vp = viewport(layout.pos.row = offset.row, layout.pos.col = 2+offset.col), hjust = -0.2, vjust = 1.7)
      grid.text(paste0(transitions[2,1],"->",transitions[2,2]), vp = viewport(layout.pos.row = offset.row, layout.pos.col = 7+offset.col), hjust = -0.2, vjust = 1.7)
      grid.text(paste0(transitions[3,1],"->",transitions[3,2]), vp = viewport(layout.pos.row = offset.row, layout.pos.col = 12+offset.col), hjust = -0.2, vjust = 1.7)
      grid.text(paste0(transitions[4,1],"->",transitions[4,2]), vp = viewport(layout.pos.row = offset.row+5, layout.pos.col = 2+offset.col), hjust = -0.2, vjust = 1.9)
      grid.text(paste0(transitions[5,1],"->",transitions[5,2]), vp = viewport(layout.pos.row = offset.row+5, layout.pos.col = 7+offset.col), hjust = -0.2, vjust = 1.9)
      grid.text(paste0(transitions[6,1],"->",transitions[6,2]), vp = viewport(layout.pos.row = offset.row+5, layout.pos.col = 12+offset.col), hjust = -0.2, vjust = 1.9)
      
      
      

      apply(transitions, 1, function(transition) {
        # transition=transitions[i,]
        # transition=transitions[6,]
        from=transition["FROM"]
        to=transition["TO"]
        fromto = paste0(from,to)
        n.total = sum(substr(data.bq.triplet$REVCOMP_TRIPLET,2,2)==as.character(from) & substr(data.bq.triplet$REVCOMP_TRIPLET,3,3)==as.character(to))
        for (baseBefore in c("A","C","G","T")) {
          for (baseAfter in c("A","C","G","T")) {
            # baseBefore="G"
            # baseAfter="G"
            triplet = paste0(baseBefore,from,to,baseAfter, collapse = "")
            transitionSubset = data.bq.triplet[data.bq.triplet$REVCOMP_TRIPLET == triplet,]
            # plot(1, xlim=c(0,45), ylim=c(0,0.4), main=paste0("Base Quality distribution for ",baseBefore,from,"->",to,baseAfter," in PID\n",PID), cex.main=0.8, xlab = "BaseQuality", ylab = "density")
            n = nrow(transitionSubset)
            if (n > 0) {
              
              # # transitionSubset.highAUC = transitionSubset[transitionSubset$AuC.alt.BQ30>0.5,]
              # VAF=transitionSubset[,"VAF"]
              # AUC=transitionSubset[,"AuC.alt.BQ30"]
              # model = lm(AUC ~ VAF)
              # # plot(model)
              # cor(VAF,AUC, method = "spearman")
              # x=cor.test(VAF,AUC)
              # # x$p.value
              
              
              color.numbers = "black"
              color.panel.border = "black"
              fontface.numbers = "plain"    
              fontsize.numbers = 1.5
              
              if (whatToPlot == "BQD_sampleIndividual") {
                
                baseScores.ref.counts = transitionSubset[,c(paste0("BQ.ref.",sapply(seq(MAX_BASE_QUALITY+1)-1, function(i) {i})))]
                baseScores.alt.counts = transitionSubset[,c(paste0("BQ.alt.",sapply(seq(MAX_BASE_QUALITY+1)-1, function(i) {i})))]
                
                baseScores.ref.counts.normalized = baseScores.ref.counts / rowSums(baseScores.ref.counts)
                baseScores.alt.counts.normalized = baseScores.alt.counts / rowSums(baseScores.alt.counts)

                baseScores.ref.counts.normalized.cumul = as.data.frame(t(apply(baseScores.ref.counts.normalized,1,cumsum)))
                baseScores.alt.counts.normalized.cumul = as.data.frame(t(apply(baseScores.alt.counts.normalized,1,cumsum)))
                baseScores.ref.counts.normalized.cumul$sample = rownames(baseScores.ref.counts.normalized.cumul)
                baseScores.alt.counts.normalized.cumul$sample = rownames(baseScores.alt.counts.normalized.cumul)
                # baseScores.ref.counts.normalized.cumul$VAF = transitionSubset$VAF
                # baseScores.alt.counts.normalized.cumul$VAF = transitionSubset$VAF

                molten.ref=melt(baseScores.ref.counts.normalized.cumul, id.vars = "sample")
                molten.ref$BQ=as.integer(gsub("BQ\\..+\\.(\\d+)$", "\\1", molten.ref$variable, perl=T))
                molten.ref$type="REF"
                molten.ref$VAF=transitionSubset[molten.ref$sample,"VAF"]
                molten.ref$nBases=transitionSubset[molten.ref$sample,"nBQ.ref"]
                molten.alt=melt(baseScores.alt.counts.normalized.cumul, id.vars = "sample")
                molten.alt$BQ=as.numeric(gsub("BQ\\..+\\.(\\d+)$", "\\1", molten.alt$variable, perl=T))
                molten.alt$type="ALT"
                molten.alt$VAF=transitionSubset[molten.alt$sample,"VAF"]
                molten.alt$nBases=transitionSubset[molten.alt$sample,"nBQ.alt"]
                molten = rbind(molten.ref, molten.alt)
                molten = molten[molten$BQ <= xlim,]
                remove(molten.ref, molten.alt)
              } else if (whatToPlot == "BQD") {
                baseScores.ref.counts = colSums(transitionSubset[,c(paste0("BQ.ref.",sapply(seq(MAX_BASE_QUALITY+1)-1, function(i) {i})))])
                baseScores.alt.counts = colSums(transitionSubset[,c(paste0("BQ.alt.",sapply(seq(MAX_BASE_QUALITY+1)-1, function(i) {i})))])
                
                baseScores.ref.counts.normalized = baseScores.ref.counts / sum(baseScores.ref.counts)
                baseScores.alt.counts.normalized = baseScores.alt.counts / sum(baseScores.alt.counts)                
                
                molten.ref=melt(baseScores.ref.counts.normalized)
                molten.ref$BQ=as.integer(gsub("BQ\\..+\\.(\\d+)$", "\\1", rownames(molten.ref), perl=T))
                molten.ref$type="REF"
                molten.alt=melt(baseScores.alt.counts.normalized)
                molten.alt$BQ=as.integer(gsub("BQ\\..+\\.(\\d+)$", "\\1", rownames(molten.alt), perl=T))
                molten.alt$type="ALT"
                molten = rbind(molten.ref, molten.alt)
                molten = molten[molten$BQ <= xlim,]
                # molten = molten[,c("sample","type","BQ","value")]

                
                AUC.BQ30.ref = round(sum(molten[with(molten, BQ <= 30 & type == 'REF'),"value"])*100,1)
                AUC.BQ30.alt = round(sum(molten[with(molten, BQ <= 30 & type == 'ALT'),"value"])*100,1)                
                
                # black | >30%:orange font | >40%:red font | >50% red font+border
                if (AUC.BQ30.alt > 30) {
                  fontface.numbers = "bold"
                  if (AUC.BQ30.alt > 40) {
                    color.numbers = "red"
                    if (AUC.BQ30.alt > 50) {
                      color.panel.border = "red"
                    }
                  } else {
                    color.numbers = "orange"
                  }
                }                   
              } else if (whatToPlot == "BQ_CoV") {
                molten.Cov.ref = melt(transitionSubset[,"CoV.ref"])
                molten.Cov.ref$type="REF"
                molten.Cov.alt = melt(transitionSubset[,"CoV.alt"])
                molten.Cov.alt$type="ALT"
                molten.CoV = rbind(molten.Cov.ref, molten.Cov.alt)
                remove(molten.Cov.ref, molten.Cov.alt)                
              }
              
              
              if (n<=3) {
                alpha = 0.2
              } else if (n<7) {
                alpha = 0.4
              } else if (n<=10) {
                alpha = 0.6
              } else {
                alpha=1.0
              }
              

              # alpha_surplusPenalty makes triplets with low overall number of base scores slightly less visible
              alpha_surplusPenalty.ref = 0.0
              if (sum(transitionSubset$nBQ.ref) < 50) {
                alpha_surplusPenalty.ref = 0.1
              }
              # lines(baseScores.ref.counts.normalized, col=add.alpha(greenColor, alpha-alpha_surplusPenalty.ref))
              
              alpha_surplusPenalty.alt = 0.0
              if (sum(transitionSubset$nBQ.alt) < 50) {
                alpha_surplusPenalty.alt = 0.1
              }
              # lines(baseScores.alt.counts.normalized, col=add.alpha(redColor, alpha-alpha_surplusPenalty.alt))
              

              
              row = as.integer(POSITIONS.BEFORE[baseBefore]+POSITIONS.FROMTO.ROW[fromto]*(4+1))
              col = as.integer(POSITIONS.AFTER[baseAfter]+POSITIONS.FROMTO.COL[fromto]*(4+1))
              vp = viewport(layout.pos.row = row, layout.pos.col = col)
              
              axis.label.size=3
              if (whatToPlot == "BQD_sampleIndividual") {

                axis.label.size.x=1
                axis.label.size.y=2
                
                
                theme = theme(legend.position="none", 
                              axis.title=element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.line=element_blank(),
                              #axis.text.x=element_text(size=axis.label.size.x), axis.text.y=element_text(size=axis.label.size.y),
                              # plot.margin = unit(c(0,-0.2,-0.35,-2.0), "cm"), 
                              plot.margin = unit(c(0,0.0,0.0,0.0), "cm"),
                              legend.margin = unit(c(0,0,0,0), "cm"),
                              panel.border = element_rect(colour = color.panel.border, fill = NA, size=0.0),
                              # panel.background = element_rect(fill = NA),
                              strip.text = element_text(size=2, lineheight=0.01),
                              panel.margin=unit(0.000, "lines"),
                              panel.grid.major = element_blank(), panel.grid.minor = element_blank())
                
                plot = ggplot(molten, aes(x=BQ,y=value, group=interaction(sample,type), colour=VAF)) + 
                  geom_line(aes(size=log2(nBases))) + scale_size(range = c(0.01, 0.1), limits = c(0.1,6)) +
                  scale_color_gradientn( colours = rainbow(7), values=c(0, 0.125, 0.25, 0.375, 0.5, 0.75, 1), limits = c(0,1), breaks=breaks) +
                  # scale_x_continuous(expand = c(-1, 0)) +
                  # scale_y_continuous(expand = c(-1, 0)) +
                  geom_errorbar(stat = "hline", yintercept = "median", width=0.005,aes(ymax=..y..,ymin=..y..), linetype=4) +
                  facet_wrap(~ type) + 
                  theme
                
                g <- ggplot_gtable(ggplot_build(plot))
                ia <- which(g$layout$name == "axis_b-1")
                g$grobs[[which(g$layout$name == "axis_b-1")]]$children[[2]] <- NULL
                g$grobs[[which(g$layout$name == "axis_b-2")]]$children[[2]] <- NULL
                g$grobs[[which(g$layout$name == "axis_l-1")]]$children[[2]] <- NULL
                g$grobs[[which(g$layout$name == "axis_l-2")]]$children[[2]] <- NULL
                grid.draw(g)
                
                g$grobs[[which(g$layout$name == "panel-1")]]$children[[5]]$children
                plot
                
                
                if (col == offset.col+1) {
                  grid.text(baseBefore, vp = viewport(layout.pos.row = row, layout.pos.col = col-1), vjust = -1.0, gp = gpar(fontsize = 6))
                }  
                if (row == offset.row+9) {
                  grid.text(baseAfter, vp = viewport(layout.pos.row = row+1, layout.pos.col = col), hjust = -0.5, gp = gpar(fontsize = 6))
                }                                 
                
                                
              }else if (whatToPlot == "BQD") {
                if (col == offset.col+1) {
                  grid.text(baseBefore, vp = viewport(layout.pos.row = row, layout.pos.col = col-1), vjust = -1.0, gp = gpar(fontsize = 6))
                  if (row == offset.row+4 | row == offset.row+9) {
                    if (row == offset.row+9) {
                      grid.text(baseAfter, vp = viewport(layout.pos.row = row+1, layout.pos.col = col), hjust = -1.0, gp = gpar(fontsize = 6))
                    }                  
                    # both axes have ticks and labels
                    theme = theme(legend.position="none", axis.title=element_blank(), 
                                  axis.ticks = element_blank(), axis.text=element_text(size=axis.label.size),
                                  plot.margin = unit(c(0,-0.1,0,-0.1), "cm"), legend.margin = unit(c(0,0,0,0), "cm"),
                                  panel.border = element_rect(colour = color.panel.border, fill = NA),
                                  panel.background = element_rect(fill = NA))
                    
                  } else {                  # only axis ticks and labels for y axis
                    
                    theme = theme(legend.position="none", axis.title=element_blank(), axis.ticks = element_blank(),
                                  plot.margin = unit(c(0,-0.05,-0.25,-0.1), "cm"), legend.margin = unit(c(0,0,0,0), "cm"),
                                  axis.text.x = element_blank(), axis.text.y=element_text(size=axis.label.size), #axis.ticks.length = unit(0,"null"), #axis.ticks.margin = unit(0,"null"),
                                  panel.border = element_rect(colour = color.panel.border, fill = NA),
                                  panel.background = element_rect(fill = NA))
                    
                  }
                } else {
                  if (row == offset.row+4 | row == offset.row+9) {
                    if (row == offset.row+9) {
                      grid.text(baseAfter, vp = viewport(layout.pos.row = row+1, layout.pos.col = col), hjust = -0.5, gp = gpar(fontsize = 6))
                    }
                    # only axis ticks and labels for x axis
                    theme = theme(legend.position="none", axis.title=element_blank(), axis.ticks = element_blank(),
                                  plot.margin = unit(c(0,-0.1,0,-0.35), "cm"), legend.margin = unit(c(0,0,0,0), "cm"),
                                  axis.text.y = element_blank(), axis.text.x=element_text(size=axis.label.size), #axis.ticks.length = unit(0,"null"), #axis.ticks.margin = unit(0,"null"),
                                  panel.border = element_rect(colour = color.panel.border, fill = NA),
                                  panel.background = element_rect(fill = NA))                  
                  } else {
                    # no axis ticks at all
                    theme = theme(legend.position="none", axis.title=element_blank(), axis.ticks = element_blank(),
                                  plot.margin = unit(c(0,-0.1,0,-0.1), "cm"), legend.margin = unit(c(0,0,0,0), "cm"),
                                  axis.text = element_blank(), axis.ticks.length = unit(0,"null"), axis.ticks.margin = unit(0,"null"),
                                  panel.border = element_rect(colour = color.panel.border, fill = NA),
                                  panel.background = element_rect(fill = NA))                  
                  }
                }                
              } else if (whatToPlot == "BQ_CoV") {
                theme = theme(legend.position="none", axis.title=element_blank(), 
                              axis.ticks = element_blank(), axis.text=element_text(size=axis.label.size),
                              plot.margin = unit(c(0,0,0,0), "cm"), legend.margin = unit(c(0,0,0,0), "cm"),
                              panel.border = element_rect(colour = color.panel.border, fill = NA),
                              panel.background = element_rect(fill = NA))
                if (col == offset.col+1) {
                  grid.text(baseBefore, vp = viewport(layout.pos.row = row, layout.pos.col = col-1), vjust = -1.0, gp = gpar(fontsize = 6))
                }  
                if (row == offset.row+9) {
                  grid.text(baseAfter, vp = viewport(layout.pos.row = row+1, layout.pos.col = col), hjust = -0.5, gp = gpar(fontsize = 6))
                }
              }

              
              if (whatToPlot == "BQD_sampleIndividual") {
                molten$sample = as.factor(molten$sample)
                molten$type = as.factor(molten$type)
                breaks <- c(0, 0.125, 0.25, 0.375, 0.5, 0.75,1)
                plot = ggplot(molten, aes(x=BQ,y=value, group=interaction(sample,type), colour=VAF)) + 
                  geom_line(aes(size=log2(nBases))) + scale_size(range = c(0.01, 0.1), limits = c(0.1,6)) +
                  scale_color_gradientn( colours = rainbow(7), values=c(0, 0.125, 0.25, 0.375, 0.5, 0.75, 1), limits = c(0,1), breaks=breaks) +
                  geom_errorbar(stat = "hline", yintercept = "median", width=0.005,aes(ymax=..y..,ymin=..y..), linetype=4) +
                  facet_wrap(~ type) + 
                  theme
              } else if (whatToPlot == "BQD") {
                
                plot = ggplot(molten, aes(x=BQ,y=value, colour=type)) + geom_line(size=0.15) + 
                  xlim(0, xlim) + ylim(0,0.4) +
                  scale_color_manual(values=c(as.character(add.alpha(redColor, alpha-alpha_surplusPenalty.alt)), as.character(add.alpha(greenColor, alpha-alpha_surplusPenalty.ref)))) +
                  theme                
              } else if (whatToPlot == "BQ_CoV") {
                molten.CoV$type = factor(molten.CoV$type, levels=c("REF","ALT"))
                plot = ggplot(molten.CoV, aes(x=type, y=value, colour=type)) + geom_boxplot(lwd=0.2, outlier.size = 0.6) + 
                  scale_y_continuous(limits=c(0,100)) +
                  scale_color_manual(values=c(as.character(add.alpha(greenColor, alpha-alpha_surplusPenalty.ref)),as.character(add.alpha(redColor, alpha-alpha_surplusPenalty.alt)))) +
                  # ylab("CoV [%]") +
                  theme(legend.position="none", axis.title.x=element_blank(), axis.title.y=element_blank(),
                        # axis.title.y=element_text(size = 3),
                        axis.ticks = element_blank(), axis.text=element_text(size=2.5),
                        plot.margin = unit(c(0,-0.1,0,-0.1), "cm"), legend.margin = unit(c(0,0,0,0), "cm"),
                        panel.border = element_rect(colour = color.panel.border, fill = NA),
                        panel.background = element_rect(fill = NA))                
                
              }

              
              # ggplot() + geom_boxplot(data=molten.CoV, aes(x=type, y=value, colour=type)) + theme

              
                                
              print(plot, vp = vp)

              if (whatToPlot == "BQD") {
                grid.text(n, vp = viewport(layout.pos.row = row, layout.pos.col = col), hjust = 0, vjust = -5.8, gp = gpar(fontsize = 2.5, col="black", fontface=fontface.numbers))
                grid.text(paste0(round(n/n.total*100,2),"%"), vp = viewport(layout.pos.row = row, layout.pos.col = col), hjust = 0.0, vjust = -4.2, gp = gpar(fontsize = 2.5, col="black", fontface=fontface.numbers))
                grid.text(paste0("ref.BQ<=30: ",AUC.BQ30.ref,"%"), vp = viewport(layout.pos.row = row, layout.pos.col = col), hjust = 0.4, vjust = -4.3, gp = gpar(fontsize = 1.50, col="black", fontface=fontface.numbers))
                grid.text(paste0("alt.BQ<=30: ",AUC.BQ30.alt,"%"), vp = viewport(layout.pos.row = row, layout.pos.col = col), hjust = 0.4, vjust = -3.0, gp = gpar(fontsize = 1.50, col=color.numbers, fontface=fontface.numbers))
                
              }

              # if (triplet == "GTGG" & n > 5) {
              #   # check if Snager artifact is detected
              #   if (AUC.BQ30.ref > 0.4 & AUC.BQ30.alt > 0.4) {
              #     d = as.data.frame(t(c("AUC.BQ30.ref"=AUC.BQ30.ref, "AUC.BQ30.alt"=AUC.BQ30.alt)))
              #     write.table(d, file = paste0(MPILEUP_FOLDER,"SangerArtifactDetected.txt"), sep = "\t", row.names = F, col.names = T, quote = F)
              #   }
              # }
            }
          }
        }
      })
      
    } else {
      if (COMBINE_REVCOMP) {
        par(mfrow=c(2,3))
      } else {
        par(mfrow=c(4,3))  
      }
      apply(transitions, 1, function(transition) {
        # transition=transitions[i,]
        # transition=transitions[1,]
        from=transition["FROM"]
        to=transition["TO"]
        
        plot(1, xlim=c(0,45), ylim=c(0,0.4), main=paste0("Base Quality distribution for ",from,"->",to," in PID\n",PID), cex.main=0.8, xlab = "BaseQuality", ylab = "density")
        
        transitionSubset = data.bq.triplet[substr(data.bq.triplet$REVCOMP_TRIPLET,2,2)==as.character(from) & substr(data.bq.triplet$REVCOMP_TRIPLET,3,3)==as.character(to),]
        if (nrow(transitionSubset) > 0) {
          for (baseBefore in c("A","C","G","T")) {
            for (baseAfter in c("A","C","G","T")) {
              # baseBefore="A"
              # baseAfter="T"
              transitionSubset.triplet = transitionSubset[transitionSubset$REVCOMP_TRIPLET==paste0(baseBefore,from,to,baseAfter),]
              n = nrow(transitionSubset.triplet)
              if (n > 0) {
                if (n<=3) {
                  alpha = 0.2
                } else if (n<7) {
                  alpha = 0.4
                } else if (n<=10) {
                  alpha = 0.6
                } else {
                  alpha=1.0
                }
                
                baseScores.ref.counts = colSums(transitionSubset.triplet[,c(paste0("BQ.ref.",sapply(seq(MAX_BASE_QUALITY+1)-1, function(i) {i})))])
                baseScores.alt.counts = colSums(transitionSubset.triplet[,c(paste0("BQ.alt.",sapply(seq(MAX_BASE_QUALITY+1)-1, function(i) {i})))])
                # baseScores.ref.counts = colSums(transitionSubset.triplet[,c(paste0("BQ.ref.",sapply(seq(50), function(i) {i})))])
                # baseScores.alt.counts = colSums(transitionSubset.triplet[,c(paste0("BQ.alt.",sapply(seq(50), function(i) {i})))])
                
                baseScores.ref.counts.normalized = baseScores.ref.counts / sum(baseScores.ref.counts)
                baseScores.alt.counts.normalized = baseScores.alt.counts / sum(baseScores.alt.counts)
                
                # alpha_surplusPenalty makes triplets with low overall number of base scores slightly less visible
                alpha_surplusPenalty.ref = 0.0
                if (sum(transitionSubset.triplet$nBQ.ref) < 50) {
                  alpha_surplusPenalty.ref = 0.1
                }
                lines(baseScores.ref.counts.normalized, col=add.alpha(greenColor, alpha-alpha_surplusPenalty.ref))
                
                alpha_surplusPenalty.alt = 0.0
                if (sum(transitionSubset.triplet$nBQ.alt) < 50) {
                  alpha_surplusPenalty.alt = 0.1
                }
                lines(baseScores.alt.counts.normalized, col=add.alpha(redColor, alpha-alpha_surplusPenalty.alt))
              }
            }
          }
          
          
          # lines(data.bq.background.dens, lwd=3)
          if (USE_BACKGROUND_BQs) {
            baseScores.ref = unlist(apply(transitionSubset, 1, function(line) {
              baseScores = as.integer(unlist(strsplit(unlist(strsplit(as.character(line["BQ_string.ref"]), ";"))[1], ",")))
            }))
            baseScores.alt = unlist(apply(transitionSubset, 1, function(line) {
              baseScores = as.integer(unlist(strsplit(unlist(strsplit(as.character(line["BQ_string.alt"]), ";"))[1], ",")))
            }))
            
            baseScores.ref.dens = density(baseScores.ref, bw = parameter.density.bandwidth, na.rm = T)
            baseScores.alt.dens = density(baseScores.alt, bw = parameter.density.bandwidth, na.rm = T)
            
            lines(data.bq.background.dens, lwd=2)
            lines(baseScores.ref.dens, col="green", lwd=2)
            lines(baseScores.alt.dens, col="red", lwd=2)
          }    
        } else {
          text(21, 0.2, "No somatic SNVS with this triplet.")
        }
      })      
    }
  dev.off()
}





if (ALT.MEDIAN.THRESHOLD > -1) {
  # PDF_OUTPUT_FILE.filtered = sub(".pdf", paste0("_filteredAltMedian",ALT.MEDIAN.THRESHOLD,".pdf"), PDF_OUTPUT_FILE)
  # plot_BQD_to_pdf(PDF_OUTPUT_FILE = PDF_OUTPUT_FILE.filtered, data.bq.triplet = data.bq.triplet.filtered)
  plot_BQD_to_pdf(PDF_OUTPUT_FILE = PDF_OUTPUT_FILE, data.bq.triplet = data.bq.triplet.filtered, whatToPlot = "BQD")
  # plot_BQD_to_pdf(PDF_OUTPUT_FILE = PDF_OUTPUT_FILE, data.bq.triplet = data.bq.triplet.filtered, whatToPlot = "BQ_CoV")

  VCF = read.table(pipe(paste0("cat ",vcfInputFile," | grep -v '^##' ")), comment.char = '', sep = "\t", header = T, stringsAsFactors = F, check.names = F)
  colnames(VCF)[1] = "CHROM"
  VCF.filtered = merge(VCF, data.bq.triplet.filtered[,c("CHROM","POS","REF","ALT")])
  chrOrder <-c((1:22),"X","Y")
  VCF.filtered$CHROM =  factor(VCF.filtered$CHROM, chrOrder, ordered=TRUE)
  VCF.filtered = VCF.filtered[with(data=VCF.filtered, order(CHROM, POS)),]
  colnames(VCF.filtered)[1] = "#CHROM"
  ncol=ncol(VCF.filtered)
  VCF.filtered = VCF.filtered[,c(1,2,5,3,4,6:ncol)]
  write.table(VCF.filtered, file=VCF_OUTPUT_FILE.filtered, quote = F, row.names = F, col.names = T, sep = "\t")  
} else {
  plot_BQD_to_pdf(PDF_OUTPUT_FILE = PDF_OUTPUT_FILE, data.bq.triplet = data.bq.triplet, whatToPlot = "BQD")
  # plot_BQD_to_pdf(PDF_OUTPUT_FILE = PDF_OUTPUT_FILE, data.bq.triplet = data.bq.triplet, whatToPlot = "BQ_CoV")
  # plot_BQD_to_pdf(PDF_OUTPUT_FILE = PDF_OUTPUT_FILE, data.bq.triplet = data.bq.triplet, whatToPlot = "BQD_sampleIndividual")
  # data.bq.triplet.short = data.bq.triplet[data.bq.triplet$REF=="T" & data.bq.triplet$ALT=="G",]
  # plot_BQD_to_pdf(PDF_OUTPUT_FILE = PDF_OUTPUT_FILE, data.bq.triplet = data.bq.triplet.short, whatToPlot = "BQD_sampleIndividual")
}




