#!/usr/bin/env Rscript
#
# Copyright (c) 2018 German Cancer Research Center (DKFZ).
#
# Distributed under the GPL-2/GPL-3 License (http://www.gnu.de/documents/gpl-2.0.en.html, https://www.gnu.org/licenses/gpl-3.0.en.html).
#

# ATTENTION: In the output files of this script multi alternative SNVs have been removed!
# This script creates QC plots
# More specific: cumulative distributions of base qualities of the somatic SNVs.
# Distribution are grouped according to the context triplet.

library(getopt)
library(ggplot2)
library(Biostrings) # for reverseComplement
library(grid) # for unit,gpar
library(reshape2)
library(RColorBrewer)


CANNEL_INDIVIDUAL_GRAPHS=T
MAX_BASE_QUALITY=55
ALT.MEDIAN.THRESHOLD = -1 # -1 means NO FILTERING
parameter.density.bandwidth = 2.00
SEQUENCE_CONTEXT_COLUMN_INDEX=11
VAF_COLUMN_INDEX=-1
SKIP_PLOTS=FALSE

opt = getopt(matrix(c(
  'vcfInputFile', 'v', 1, "character",
  'mpileupFolder', 'm', 1, "character",
  'alignmentFolder', 'a', 1, "character",
  'PID', 'p', 1, "character",
  'background', 'b', 2, "integer",
  'outFilePrefix', 'o', 1, "character",
  'forceRerun', 'R', 2, "integer",
  'combineRevcomp', 'c', 2, "integer",
  'filterThreshold', 'f', 2, "integer",
  'sequenceContextColumnIndex', 's', 2, "integer",
  'MAFColumnIndex', 'z', 1, "integer",
  'channelIndividualGraphs', 'i', 2, "integer",
  'mainTitle', 't', 2, "character",
  'skipPlots', 'x', 2, "integer",
  'refBaseQual', 'k', 1, "character",
  'altBaseQual', 'l', 1, "character",
  'altReadPos', 'q', 1, "character",
  'refReadPos', 'r', 1, "character"
),ncol=4,byrow=TRUE));



checkForMissingParameter = function(parameter, errorText, exitCode=1) {
  if (is.null(opt[[parameter]])){
    cat(paste0(errorText,"\n")) 
    q(save = "no", status = exitCode)      # quit, status unequal 0 means error
  }
}

setIntegerValueFromParameter = function(parameterName, valueName) {
  if (! is.null(opt[[parameterName]])){
    tmp = as.integer(opt[[parameterName]])
    if (!is.na(tmp)) {
      assign(valueName, tmp, envir = .GlobalEnv)
    }
  }  
}


checkForMissingParameter("vcfInputFile", "Please specify the file that contains the SNVs for which the base score distribution plot shall be created.", 1)
    vcfInputFile = paste0(opt$vcfInputFile)
checkForMissingParameter("mpileupFolder", "Please specify the mpileup folder.", 1)
    MPILEUP_FOLDER=paste0(opt$mpileupFolder,"/")
checkForMissingParameter("alignmentFolder", "Please specify the alignment folder.", 1)
    ALIGNMENT_FOLDER=paste0(opt$alignmentFolder,"/")
checkForMissingParameter("PID", "Please specify the PID.", 1)
    PID=opt$PID
checkForMissingParameter("refBaseQual", "Please specify the reference allel base qualities file.", 1)    
    RefAlleleBaseQualitiesFile=opt$refBaseQual
checkForMissingParameter("altBaseQual", "Please specify the alternative allel base qualities file.", 1)    
    AltAlleleBaseQualitiesFile=opt$altBaseQual
checkForMissingParameter("altReadPos", "Please specify the alternative allele read positions file.", 1)    
    AltAlleleReadPositionsFile=opt$altReadPos
checkForMissingParameter("refReadPos", "Please specify the reference allele read positions file.", 1)    
    RefAlleleReadPositionsFile=opt$refReadPos
checkForMissingParameter("outFilePrefix", "Please specify the output pdf file.", 1)
    PDF_OUTPUT_FILE_PREFIX=paste0(opt$outFilePrefix)


if (is.null(opt$background)){
  USE_BACKGROUND_BQs = FALSE
} else {
  USE_BACKGROUND_BQs = ifelse(opt$background==1, TRUE, FALSE)
}

if (is.null(opt$forceRerun)){     
  FORCE_RERUN = FALSE
} else {
  FORCE_RERUN = ifelse(opt$forceRerun==1, TRUE, FALSE)
}

if (is.null(opt$combineRevcomp)){
  COMBINE_REVCOMP = TRUE
} else {
  COMBINE_REVCOMP = ifelse(opt$combineRevcomp==1, TRUE, FALSE)
}

if (! is.null(opt$channelIndividualGraphs)){
  CANNEL_INDIVIDUAL_GRAPHS = ifelse (opt$channelIndividualGraphs == 1, TRUE, FALSE)
}

if (! is.null(opt$skipPlots)){
  SKIP_PLOTS = ifelse (opt$skipPlots == 1, TRUE, FALSE)
}


setIntegerValueFromParameter(parameterName = "MAFColumnIndex", valueName = "VAF_COLUMN_INDEX")
setIntegerValueFromParameter(parameterName = "filterThreshold", valueName = "ALT.MEDIAN.THRESHOLD")
setIntegerValueFromParameter(parameterName = "sequenceContextColumnIndex", valueName = "SEQUENCE_CONTEXT_COLUMN_INDEX")

# Define possible transitions
if (COMBINE_REVCOMP) {
  # combine SNVs of a certain triplet class with the SNVs of the reverse complement triplet class?
  transitions = data.frame( c(rep("C",3),rep("T",3)), c("A","G","T","A","C","G"), stringsAsFactors = F)  
  possible_mutations = c("CA", "CG", "CT", "TA", "TC", "TG")
  forbidden_mutations = sapply(possible_mutations, function(mutation) {as.character(complement(DNAString(mutation)))})
} else {
  transitions = data.frame( c(rep("A",3),rep("C",3),rep("G",3),rep("T",3)), c("C","G","T","A","G","T","A","C","T","A","C","G"), stringsAsFactors = F)
  possible_mutations = c("AC","AG","AT", "CA","CG","CT", "GA","GC","GT", "TA","TC","TG")
}
colnames(transitions) = c("FROM", "TO")




# ALT.MEDIAN.THRESHOLD
if (ALT.MEDIAN.THRESHOLD > -1) {
  DATA_RESULTS_FILE=paste0(MPILEUP_FOLDER,"BQ_TripletSpecific_MedianFiltered",ALT.MEDIAN.THRESHOLD,".RData")
  VCF_OUTPUT_FILE.filtered = sub(".vcf$", paste0("_filteredAltMedian",ALT.MEDIAN.THRESHOLD,".vcf"), vcfInputFile)  
} else {
  DATA_RESULTS_FILE=paste0(MPILEUP_FOLDER,"BQ_TripletSpecific.RData")
}


if ( file.exists(DATA_RESULTS_FILE) & ! FORCE_RERUN ) {
  cat (paste0("Loading triplet-specific data from file ",DATA_RESULTS_FILE,"\n"))
  lnames = load(file = DATA_RESULTS_FILE)
  if ( USE_BACKGROUND_BQs & !exists("data.bq.background.dens") ) {
    data.bq.background = read.table(paste0(ALIGNMENT_FOLDER,"BaseQualitiesForChosenPositions.1000000.1.txt.gz", collapse = ""))
    colnames(data.bq.background) = "BQ"
    data.bq.background.dens = density(data.bq.background$BQ, bw = parameter.density.bandwidth)
    remove(data.bq.background)
  }
} else {
  cat (paste0("Extracting triplet-specific data...\n"))
  # GO THROUGH DATA PROCESSING STEPS (either due to missing data file (e.g. first run) or being forced to do it although data file is available)
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
  cat (paste0("Finished reading vcf file...\n"))
  data.bq.triplet$Triplet = with(data.bq.triplet, paste0(BaseBefore,REF,ALT,BaseAfter))

  if (length(grep(',', data.bq.triplet$"Triplet"))>0) {
    data.bq.triplet = data.bq.triplet[-grep(',', data.bq.triplet$"Triplet"),]
  }
  
  data.bq.ref = read.table(file = RefAlleleBaseQualitiesFile, stringsAsFactors = F)
  data.bq.alt = read.table(file = AltAlleleBaseQualitiesFile, stringsAsFactors = F)
  data.rp.alt = read.table(file = AltAlleleReadPositionsFile, stringsAsFactors = F)
  
  colnames(data.bq.ref) = c("CHROM","POS","BQ_string")
  colnames(data.bq.alt) = c("CHROM","POS","BQ_string")
  colnames(data.rp.alt) = c("CHROM","POS","RP_string")
  
  cat ("Merging Base Qualities and Read Positions to dataframe...\n")
  data.bq.triplet = merge(data.bq.triplet, data.bq.ref, by=c("CHROM","POS"), all.x = T)
  colnames(data.bq.triplet)[ncol(data.bq.triplet)] = "BQ_string.ref"
  data.bq.triplet = merge(data.bq.triplet, data.bq.alt, by=c("CHROM","POS"), all.x = T)
  colnames(data.bq.triplet)[ncol(data.bq.triplet)] = "BQ_string.alt"
  data.bq.triplet = merge(data.bq.triplet, data.rp.alt, by=c("CHROM","POS"), all.x = T)
  colnames(data.bq.triplet)[ncol(data.bq.triplet)] = "RP_string.alt"
  data.bq.triplet = data.bq.triplet[order(data.bq.triplet$Triplet),]
  
  cat (paste0("Generating triplet-specific base quality matrix...\n"))
  data.bq.triplet.processed = lapply(seq(nrow(data.bq.triplet)), function(i) {
    line = data.bq.triplet[i,]
    baseScores.ref = as.integer(unlist(strsplit(unlist(strsplit(as.character(line[1,"BQ_string.ref"]), ";"))[1], ",")))
    baseScores.alt = as.integer(unlist(strsplit(unlist(strsplit(as.character(line[1,"BQ_string.alt"]), ";"))[1], ",")))
    readPos.alt = as.integer(unlist(strsplit(unlist(strsplit(as.character(line[1,"RP_string.alt"]), ";"))[1], ",")))
    
    countsLine.ref = data.frame(matrix(0, ncol = MAX_BASE_QUALITY+1, nrow = 1))
    colnames(countsLine.ref) = c(paste0("BQ.ref.",sapply(seq(MAX_BASE_QUALITY+1)-1, function(i) {i})))
    countsLine.alt = data.frame(matrix(0, ncol = MAX_BASE_QUALITY+1, nrow = 1))
    colnames(countsLine.alt) = c(paste0("BQ.alt.",sapply(seq(MAX_BASE_QUALITY+1)-1, function(i) {i})))    
    counts.ref = table(baseScores.ref);  if (length(counts.ref)>0) {countsLine.ref[1,paste0("BQ.ref.",names(counts.ref))] = as.integer(unlist(counts.ref))}
    if(length(counts.ref)==0) {counts.ref=NA}
    counts.alt = table(baseScores.alt);  if (length(counts.alt)>0) {countsLine.alt[1,paste0("BQ.alt.",names(counts.alt))] = as.integer(unlist(counts.alt))}
    if(length(counts.alt)==0) {counts.alt=NA}
    
    line[1,"nBQ.ref"] = ifelse(length(baseScores.ref)==1 && is.na(baseScores.ref),0,length(baseScores.ref))
    line[1,"nBQ.alt"] = ifelse(length(baseScores.alt)==1 && is.na(baseScores.alt),0,length(baseScores.alt))
    
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
  cat (paste0("Finished generation of triplet-specific base quality matrix...\n"))
  
  data.bq.triplet$"BaseBefore" = factor(data.bq.triplet$"BaseBefore", levels=c("A","C","G","T"))
  data.bq.triplet$"REF" = factor(data.bq.triplet$"REF", levels=c("A","C","G","T"))
  data.bq.triplet$"ALT" = factor(data.bq.triplet$"ALT", levels=c("A","C","G","T"))
  data.bq.triplet$"BaseAfter" = factor(data.bq.triplet$"BaseAfter", levels=c("A","C","G","T"))
  if ( ! is.null(data.bq.triplet$VAF) ) {
    data.bq.triplet$"VAF" = as.numeric(as.character(data.bq.triplet$"VAF"))
  }
  data.bq.triplet$"Triplet" = as.factor(as.character(data.bq.triplet$"Triplet"))
  
  
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
  
  
  if (USE_BACKGROUND_BQs) {
    save(data.bq.triplet, data.bq.background.dens, file = DATA_RESULTS_FILE)
    # write.table(data.bq.triplet.filtered)
  } else {
    save(data.bq.triplet, file = DATA_RESULTS_FILE)
  }
  cat (paste0("Finished extraction of triplet-specific data...\n"))
}



# Filter SNVs by median ALT BQ
# data.bq.triplet.filtered will contain only those SNVs whose median of alternative-allele base qualities is at least ALT.MEDIAN.THRESHOLD
if (ALT.MEDIAN.THRESHOLD > -1) {
  data.bq.triplet.filtered = do.call(rbind, apply(transitions, 1, function(TRANSITION) {
    if ( "THIS SNV IS ARTIFACTED" == "THIS SNV IS ARTIFACTED" ) {
      # think about: only filter channels that show artifact?
      SNVs = data.bq.triplet[with(data.bq.triplet, substr(REVCOMP_TRIPLET,2,2)==TRANSITION[1] & substr(REVCOMP_TRIPLET,3,3)==TRANSITION[2]),]
      indices = which(SNVs$MeadianBQ.alt >= ALT.MEDIAN.THRESHOLD)
      # indices = which(SNVs$MeadianBQ.alt >= ALT.MEDIAN.THRESHOLD | SNVs$nBQ.alt <= 5)
      # SNVs$nBQ.alt
      SNVs = SNVs[indices,]    
      return (SNVs)
    }
  }))
}


######################################################
# PLOTTING PART BEGINS HERE
######################################################

cat (paste0("Starting plot generation...\n"))

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






scatterplot_altBQ_vs_altReadPos_forTriplet = function(data.bq.triplet, triplet, AUC_threshold=AUC_threshold) {
# plot base quality (alternative allele) vs. read position for a single triplet
  subset = data.bq.triplet[data.bq.triplet$REVCOMP_TRIPLET == triplet,]
  plot(-5, xlab="Read Position", ylab="Base Quality", xlim=c(0,110), ylim=c(1,50), main=triplet)
  apply(subset, 1, function(line) {
    BQ = as.integer(unlist(strsplit(unlist(strsplit(as.character(line["BQ_string.alt"]),";"))[1],",")))
    ReadPos = as.integer(unlist(strsplit(unlist(strsplit(as.character(line["RP_string.alt"]),";"))[1],",")))
    AUC.alt.BQ30 = as.numeric(line["AuC.alt.BQ30"])
    color = ifelse(AUC.alt.BQ30 >= AUC_threshold, "red", "black")
    pch = ifelse(AUC.alt.BQ30 >= AUC_threshold, 20, ".")
    points(ReadPos, BQ, pch = pch, col=color)
  })
}

plot_altBQ_vs_altReadPos = function(data.bq.triplet, transitions, AUC_threshold=0.25, RESULT_PDF = NA) {
# plot base quality (alternative allele) vs. read position for all triplets (calls single triplet plotter function)
  if (!is.na(RESULT_PDF)) {pdf(RESULT_PDF)}
  par(mfrow=c(4,4))
  for (i in seq(nrow(transitions))) {
    transition = transitions[i,]
    for (baseBefore in c("A","C","G","T")) {
      for (baseAfter in c("A","C","G","T")) {
        triplet = paste0(baseBefore,transition["FROM"],transition["TO"],baseAfter, collapse = "")
        scatterplot_altBQ_vs_altReadPos_forTriplet(data.bq.triplet = data.bq.triplet, triplet = triplet, AUC_threshold = AUC_threshold)
      }
    }
  }
  if (!is.na(RESULT_PDF)) {dev.off()}
}


plot_BQD_to_pdf = function(PDF_OUTPUT_FILE, data.bq.triplet, 
                           whatToPlot=c("BQD","BQ_CoV","BQD_sampleIndividual","BQD_sampleIndividual_ColoredByChromosome", "BQD_sampleIndividual_ColoredByVAF","BQD_sampleIndividual_ColoredByReadPosition"), 
                           ReadPositionQuantile=0.75) {
# function to plot cumulative base quality distributions with different color schemes
# triplets with a AUC.BQ30.alt (area under the curve BQ 30 = proportion of alternative allele reads with base quality <= 30) value above 30%/40%/50% will be highlighted in the plots
  if (! is.na(PDF_OUTPUT_FILE)) {
    print(paste0("Creating output pdf file ",PDF_OUTPUT_FILE))
    pdf(PDF_OUTPUT_FILE, 8.27, 6.0)
    # pdf(PDF_OUTPUT_FILE, paper = "a4r")
    # par(mar=c(1,1,1,1)+0.1, mar=c(1,1,1,1))
  }
  ReadPositionsQuantile_ColName = paste0("ReadPositions.Q",ReadPositionQuantile*100)
  
  
    if (CANNEL_INDIVIDUAL_GRAPHS == T) {
      if (COMBINE_REVCOMP == F) {
        stop("Channel-individual plots are only supported for combined revcomp!")
      }
      
      grid.newpage(recording = F)
      if (whatToPlot == "BQD") {
        gl <- grid.layout(nrow=13, ncol=17, widths = unit(c(3,1,2.6,2,2,2,1,2,2,2,2,1,2,2,2,2,3), "null"), heights = unit(c(2,2,2,2,2,2.5,2.5,2,2,2,2.5,0.5), "null"))
      } else if (whatToPlot == "BQ_CoV") {
        gl <- grid.layout(nrow=13, ncol=17, widths = unit(c(3,1,2,2,2,2,1,2,2,2,2,1,2,2,2,2,3), "null"), heights = unit(c(2,2,2,2,2,2,2.5,2,2,2,2,0.5), "null"))
      } else if (whatToPlot == "BQD_sampleIndividual" || whatToPlot == "BQD_sampleIndividual_ColoredByChromosome" || whatToPlot == "BQD_sampleIndividual_ColoredByVAF" || whatToPlot == "BQD_sampleIndividual_ColoredByReadPosition") {
        gl <- grid.layout(nrow=13, ncol=19, widths = unit(c(2,1,2,2,2,2,1,2,2,2,2,1,2,2,2,2,2,2,1), "null"), heights = unit(c(2,2,2,2,2,2,2.5,2,2,2,2,0.5), "null"))
        legend.x=18
        legend.y=c(5,9)
      }
      
      pushViewport(viewport(layout = gl))
      
      grid.text(opt$mainTitle, vp = viewport(layout.pos.row = 1, layout.pos.col = 2+offset.col), hjust = -0.2, vjust = 1.7)
      
      grid.text(paste0(transitions[1,1],"->",transitions[1,2]), vp = viewport(layout.pos.row = offset.row, layout.pos.col = 2+offset.col), hjust = -0.2, vjust = 1.7)
      grid.text(paste0(transitions[2,1],"->",transitions[2,2]), vp = viewport(layout.pos.row = offset.row, layout.pos.col = 7+offset.col), hjust = -0.2, vjust = 1.7)
      grid.text(paste0(transitions[3,1],"->",transitions[3,2]), vp = viewport(layout.pos.row = offset.row, layout.pos.col = 12+offset.col), hjust = -0.2, vjust = 1.7)
      grid.text(paste0(transitions[4,1],"->",transitions[4,2]), vp = viewport(layout.pos.row = offset.row+5, layout.pos.col = 2+offset.col), hjust = -0.2, vjust = 1.9)
      grid.text(paste0(transitions[5,1],"->",transitions[5,2]), vp = viewport(layout.pos.row = offset.row+5, layout.pos.col = 7+offset.col), hjust = -0.2, vjust = 1.9)
      grid.text(paste0(transitions[6,1],"->",transitions[6,2]), vp = viewport(layout.pos.row = offset.row+5, layout.pos.col = 12+offset.col), hjust = -0.2, vjust = 1.9)
      
      
      
      legend = NA
      apply(transitions, 1, function(transition) {
        # transition = transitions[1,]
        from=transition["FROM"]
        to=transition["TO"]
        fromto = paste0(from,to)
        n.total = sum(substr(data.bq.triplet$REVCOMP_TRIPLET,2,2)==as.character(from) & substr(data.bq.triplet$REVCOMP_TRIPLET,3,3)==as.character(to))
        for (baseBefore in c("A","C","G","T")) {
          for (baseAfter in c("A","C","G","T")) {
            triplet = paste0(baseBefore,from,to,baseAfter, collapse = "")
            transitionSubset = data.bq.triplet[data.bq.triplet$REVCOMP_TRIPLET == triplet,]
            n = nrow(transitionSubset)
            if (n > 0) {
              
              color.numbers = "black"
              color.panel.border = "black"
              fontface.numbers = "plain"    
              fontsize.numbers = 1.5
              
              if (whatToPlot == "BQD_sampleIndividual" || whatToPlot == "BQD_sampleIndividual_ColoredByChromosome" || 
                  whatToPlot == "BQD_sampleIndividual_ColoredByVAF" || whatToPlot == "BQD_sampleIndividual_ColoredByReadPosition") {
                
                baseScores.ref.counts = transitionSubset[,c(paste0("BQ.ref.",sapply(seq(MAX_BASE_QUALITY+1)-1, function(i) {i})))]
                baseScores.alt.counts = transitionSubset[,c(paste0("BQ.alt.",sapply(seq(MAX_BASE_QUALITY+1)-1, function(i) {i})))]
                
                baseScores.ref.counts.normalized = baseScores.ref.counts / rowSums(baseScores.ref.counts)
                baseScores.alt.counts.normalized = baseScores.alt.counts / rowSums(baseScores.alt.counts)

                baseScores.ref.counts.normalized.cumul = as.data.frame(t(apply(baseScores.ref.counts.normalized,1,cumsum)))
                baseScores.alt.counts.normalized.cumul = as.data.frame(t(apply(baseScores.alt.counts.normalized,1,cumsum)))
                baseScores.ref.counts.normalized.cumul$sample = rownames(baseScores.ref.counts.normalized.cumul)
                baseScores.alt.counts.normalized.cumul$sample = rownames(baseScores.alt.counts.normalized.cumul)

                molten.ref=melt(baseScores.ref.counts.normalized.cumul, id.vars = "sample")
                molten.ref$BQ=as.integer(gsub("BQ\\..+\\.(\\d+)$", "\\1", molten.ref$variable, perl=T))
                molten.ref$type="REF"
                molten.ref$VAF=transitionSubset[molten.ref$sample,"VAF"]
                molten.ref$AUC30=transitionSubset[molten.ref$sample,"AuC.ref.BQ30"]
                molten.ref$nBases=transitionSubset[molten.ref$sample,"nBQ.ref"]
                molten.ref[,ReadPositionsQuantile_ColName] = -1
                
                molten.alt=melt(baseScores.alt.counts.normalized.cumul, id.vars = "sample")
                molten.alt$BQ=as.numeric(gsub("BQ\\..+\\.(\\d+)$", "\\1", molten.alt$variable, perl=T))
                molten.alt$type="ALT"
                molten.alt$VAF=transitionSubset[molten.alt$sample,"VAF"]
                molten.alt$AUC30=transitionSubset[molten.alt$sample,"AuC.alt.BQ30"]
                molten.alt$nBases=transitionSubset[molten.alt$sample,"nBQ.alt"]
                molten.alt[,ReadPositionsQuantile_ColName] = sapply(transitionSubset[molten.alt$sample,"RP_string.alt"], function(positionsString) {
                  positions = as.integer(unlist(strsplit(positionsString, ",")))
                  return (round(quantile(positions, ReadPositionQuantile)))
                })                
                # molten.altReadPos = melt(baseScores.alt.counts.normalized.cumul, id.vars = "sample")
                if  (whatToPlot == "BQD_sampleIndividual_ColoredByChromosome") {
                  molten.ref$CHROM = factor(transitionSubset[molten.ref$sample,"CHROM"], levels = c(seq(22),"X","Y"))
                  molten.alt$CHROM = factor(transitionSubset[molten.alt$sample,"CHROM"], levels = c(seq(22),"X","Y"))
                }                
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
                
                # Define highlighting of triplets due to high number of reads with low quality alternative base calls
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
              
              # Make density plot lines less visible for a certain triplet if only few SNVs represent this triplet
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
              
              alpha_surplusPenalty.alt = 0.0
              if (sum(transitionSubset$nBQ.alt) < 50) {
                alpha_surplusPenalty.alt = 0.1
              }
              
              row = as.integer(POSITIONS.BEFORE[baseBefore]+POSITIONS.FROMTO.ROW[fromto]*(4+1))
              col = as.integer(POSITIONS.AFTER[baseAfter]+POSITIONS.FROMTO.COL[fromto]*(4+1))
              vp = viewport(layout.pos.row = row, layout.pos.col = col)
              
              axis.label.size=3
              if (whatToPlot == "BQD_sampleIndividual" || whatToPlot == "BQD_sampleIndividual_ColoredByChromosome" || 
                  whatToPlot == "BQD_sampleIndividual_ColoredByVAF" || whatToPlot == "BQD_sampleIndividual_ColoredByReadPosition") {

                axis.label.size.x=1.3
                axis.label.size.y=1.7
                
                theme = theme(legend.position="none",
                              axis.title=element_blank(),
                              axis.ticks=element_blank(),
                              axis.text.x=element_text(size=axis.label.size.x), axis.text.y=element_text(size=axis.label.size.y),
                              axis.ticks.margin = unit(-0.14, "cm"),
                              plot.margin = unit(c(0,0,0,0), "cm"),
                              panel.border = element_rect(colour = color.panel.border, fill = NA, size=0.0),
                              panel.background = element_rect(fill = NA),
                              strip.text = element_text(size=1.5, lineheight=0.15, vjust = 0.01),
                              strip.background=element_blank(),
                              panel.margin=unit(0.000, "lines"),
                              panel.grid.major = element_blank(), panel.grid.minor = element_blank())
                

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

              
              if (whatToPlot == "BQD_sampleIndividual" || whatToPlot == "BQD_sampleIndividual_ColoredByChromosome" || 
                  whatToPlot == "BQD_sampleIndividual_ColoredByVAF" || whatToPlot == "BQD_sampleIndividual_ColoredByReadPosition") {
                molten$sample = as.factor(molten$sample)
                molten$type = factor(molten$type, levels=c("REF","ALT"))
                molten$log2_nBases = log2(molten$nBases)
                  
                breaks <- c(0, 0.125, 0.25, 0.375, 0.5, 0.75,1)
                plot = ggplot(molten, aes(x=BQ,y=value, group=interaction(sample,type)))
                if (whatToPlot == "BQD_sampleIndividual_ColoredByChromosome") {
                  plot = plot + geom_line(aes(size=log2(nBases), colour=CHROM)) + scale_size(range = c(0.01, 0.1), limits = c(0.1,15)) +
                    scale_color_manual(values=rainbow(24), limits = c(seq(22),"X","Y"))
                } else if (whatToPlot == "BQD_sampleIndividual_ColoredByVAF") {
                  if (VAF_COLUMN_INDEX > -1) {
                    plot = plot + geom_line(aes(size=log2(nBases), colour=VAF)) + scale_size(range = c(0.01, 0.1), limits = c(0.1,15)) +
                      scale_color_gradientn( colours = rev(rainbow(7)), limits = c(0,1), values=breaks, breaks=breaks)
                  } else {
                    plot = plot + geom_line(aes(size=log2(nBases), colour=type)) + scale_size(range = c(0.01, 0.1), limits = c(0.1,15)) + #scale_size(range = c(0.01, 0.1)
                      scale_color_manual(values=c(greenColor,redColor))
                  }                  
                } else if (whatToPlot == "BQD_sampleIndividual") {
                  plot = plot + geom_line(aes(size=log2(nBases), colour=type)) + scale_size(range = c(0.01, 0.1), limits = c(0.1,15)) + #scale_size(range = c(0.01, 0.1)
                    scale_color_manual(values=c(greenColor,redColor))                  
                } else if (whatToPlot == "BQD_sampleIndividual_ColoredByReadPosition") {
                  plot = plot + geom_line(aes_string(size="log2_nBases", colour=ReadPositionsQuantile_ColName)) + scale_size(range = c(0.01, 0.1), limits = c(0.1,15)) +
                    scale_color_gradientn( colours = rev(rainbow(7)), limits = c(0,max(molten[,ReadPositionsQuantile_ColName])))
                }
                plot = plot + 
                  #geom_hline(aes(yintercept=AUC30), size=0.005, linetype=1) +
                  geom_vline(xintercept=30, size=0.02, linetype=1) +
                  facet_wrap(~ type) +
                  theme(legend.text=element_text(size=6), legend.key.height=unit(10,"pt"), legend.title=element_text(size=5))
                  
                #Extract Legend 
                if (!"grob" %in% class(legend)) {
                  tmp <- ggplot_gtable(ggplot_build(plot)) 
                  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
                  legend <<- tmp$grobs[[leg]]
                  pushViewport(viewport(layout.pos.row = legend.y, layout.pos.col = legend.x))
                  grid.draw(legend)
                  popViewport()
                }
                
                plot = plot + theme
                g <- ggplot_gtable(ggplot_build(plot))
                # changes necessary due to changes in ggplot2 as used in Rv3.3.1 ?
                if (as.integer(R.Version()["major"]) >= 3 & as.integer(unlist(strsplit(as.character(R.Version()["minor"]),"\\."))[1]) >= 3) {
                  t = subset(g$layout, name == "strip-t-1-1")$t
                  l = (subset(g$layout, name == "axis-l-1-1")$l)-1
                  b = subset(g$layout, name == "axis-b-1-1")$b
                  r = subset(g$layout, name == "panel-2-1")$r
                } else {
                  t = subset(g$layout, name == "strip_t-1")$t
                  l = (subset(g$layout, name == "axis_l-1")$l)-1
                  b = subset(g$layout, name == "axis_b-1")$b
                  r = subset(g$layout, name == "panel-2")$r
                }
                panel  = g[t:b, l:r]
                panel$vp = vp
                grid.draw(panel)                
              } else if (whatToPlot == "BQD") {
                
                plot = ggplot(molten, aes(x=BQ,y=value, colour=type)) + geom_line(size=0.15) + 
                  xlim(0, xlim) + ylim(0,0.4) +
                  scale_color_manual(values=c(as.character(add.alpha(redColor, alpha-alpha_surplusPenalty.alt)), as.character(add.alpha(greenColor, alpha-alpha_surplusPenalty.ref)))) +
                  theme   
                print(plot, vp = vp)
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
                print(plot, vp = vp)
              }

              
              if (whatToPlot == "BQD") {
                grid.text(n, vp = viewport(layout.pos.row = row, layout.pos.col = col), hjust = 0, vjust = -5.8, gp = gpar(fontsize = 2.5, col="black", fontface=fontface.numbers))
                grid.text(paste0(round(n/n.total*100,2),"%"), vp = viewport(layout.pos.row = row, layout.pos.col = col), hjust = 0.0, vjust = -4.2, gp = gpar(fontsize = 2.5, col="black", fontface=fontface.numbers))
                grid.text(paste0("ref.BQ<=30: ",AUC.BQ30.ref,"%"), vp = viewport(layout.pos.row = row, layout.pos.col = col), hjust = 0.4, vjust = -4.3, gp = gpar(fontsize = 1.50, col="black", fontface=fontface.numbers))
                grid.text(paste0("alt.BQ<=30: ",AUC.BQ30.alt,"%"), vp = viewport(layout.pos.row = row, layout.pos.col = col), hjust = 0.4, vjust = -3.0, gp = gpar(fontsize = 1.50, col=color.numbers, fontface=fontface.numbers))
                
              }
            }
          }
        }
        ################### 
      }) # end of apply
      
    } else {
      if (COMBINE_REVCOMP) {
        par(mfrow=c(2,3))
      } else {
        par(mfrow=c(4,3))  
      }
      apply(transitions, 1, function(transition) {
        from=transition["FROM"]
        to=transition["TO"]
        
        plot(1, xlim=c(0,45), ylim=c(0,0.4), main=paste0("Base Quality distribution for ",from,"->",to," in PID\n",PID), cex.main=0.8, xlab = "BaseQuality", ylab = "density")
        
        transitionSubset = data.bq.triplet[substr(data.bq.triplet$REVCOMP_TRIPLET,2,2)==as.character(from) & substr(data.bq.triplet$REVCOMP_TRIPLET,3,3)==as.character(to),]
        if (nrow(transitionSubset) > 0) {
          for (baseBefore in c("A","C","G","T")) {
            for (baseAfter in c("A","C","G","T")) {
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

                baseScores.ref.counts.normalized = baseScores.ref.counts / sum(baseScores.ref.counts)
                baseScores.alt.counts.normalized = baseScores.alt.counts / sum(baseScores.alt.counts)
                
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
  if (! SKIP_PLOTS) {
    print("Plotting meadian-filtered data")
    plot_BQD_to_pdf(PDF_OUTPUT_FILE = paste0(PDF_OUTPUT_FILE_PREFIX,".BQD_combined.pdf"), data.bq.triplet = data.bq.triplet.filtered, whatToPlot = "BQD")
    plot_BQD_to_pdf(PDF_OUTPUT_FILE = paste0(PDF_OUTPUT_FILE_PREFIX,".BQD_CoV.pdf"), data.bq.triplet = data.bq.triplet.filtered, whatToPlot = "BQ_CoV")
    plot_BQD_to_pdf(PDF_OUTPUT_FILE = paste0(PDF_OUTPUT_FILE_PREFIX,".BQD_individual.pdf"), data.bq.triplet = data.bq.triplet.filtered, whatToPlot = "BQD_sampleIndividual")
    plot_BQD_to_pdf(PDF_OUTPUT_FILE = paste0(PDF_OUTPUT_FILE_PREFIX,".BQD_individual_CHROMcolored.pdf"), data.bq.triplet = data.bq.triplet.filtered, whatToPlot = "BQD_sampleIndividual_ColoredByChromosome")
    plot_BQD_to_pdf(PDF_OUTPUT_FILE = paste0(PDF_OUTPUT_FILE_PREFIX,".BQD_individual_VAFcolored.pdf"), data.bq.triplet = data.bq.triplet.filtered, whatToPlot = "BQD_sampleIndividual_ColoredByVAF")

    ReadPositionQuantile=0.5
    plot_BQD_to_pdf(PDF_OUTPUT_FILE = paste0(PDF_OUTPUT_FILE_PREFIX,".BQD_individual_ReadPosColored_Q",ReadPositionQuantile*100,".pdf"), data.bq.triplet = data.bq.triplet.filtered,
    whatToPlot = "BQD_sampleIndividual_ColoredByReadPosition", ReadPositionQuantile=ReadPositionQuantile )

    ReadPositionQuantile=0.6
    plot_BQD_to_pdf(PDF_OUTPUT_FILE = paste0(PDF_OUTPUT_FILE_PREFIX,".BQD_individual_ReadPosColored_Q",ReadPositionQuantile*100,".pdf"), data.bq.triplet = data.bq.triplet.filtered,
    whatToPlot = "BQD_sampleIndividual_ColoredByReadPosition", ReadPositionQuantile=ReadPositionQuantile )

    ReadPositionQuantile=0.7
    plot_BQD_to_pdf(PDF_OUTPUT_FILE = paste0(PDF_OUTPUT_FILE_PREFIX,".BQD_individual_ReadPosColored_Q",ReadPositionQuantile*100,".pdf"), data.bq.triplet = data.bq.triplet.filtered,
    whatToPlot = "BQD_sampleIndividual_ColoredByReadPosition", ReadPositionQuantile=ReadPositionQuantile )

    ReadPositionQuantile=0.8
    plot_BQD_to_pdf(PDF_OUTPUT_FILE = paste0(PDF_OUTPUT_FILE_PREFIX,".BQD_individual_ReadPosColored_Q",ReadPositionQuantile*100,".pdf"), data.bq.triplet = data.bq.triplet.filtered,
    whatToPlot = "BQD_sampleIndividual_ColoredByReadPosition", ReadPositionQuantile=ReadPositionQuantile )

    print("BQ vs. ReadPosition Scatterplot")
    plot_altBQ_vs_altReadPos(data.bq.triplet = data.bq.triplet.filtered,
    transitions = transitions, AUC_threshold=0.25,
    RESULT_PDF = paste0(MPILEUP_FOLDER,"Scatterplot_BQ_ReadPos.pdf"))
  }

# write out the filtered file named ${SNV_FILE_WITH_MAF_filtered}
  VCF = read.table(pipe(paste0("cat ",vcfInputFile," | grep -v '^##' ")), comment.char = '', sep = "\t", header = T, stringsAsFactors = F, check.names = F)
  colnames(VCF)[1] = "CHROM"
  VCF.filtered = merge(VCF, data.bq.triplet.filtered[,c("CHROM","POS","REF","ALT")])
  chrOrder <-c((1:22),"X","Y")
  VCF.filtered$CHROM =  factor(VCF.filtered$CHROM, chrOrder, ordered=TRUE)
  VCF.filtered = VCF.filtered[with(data=VCF.filtered, order(CHROM, POS)),]
  colnames(VCF.filtered)[1] = "#CHROM"
  ncol=ncol(VCF.filtered)
  VCF.filtered = VCF.filtered[,c(1,2,5,3,4,6:ncol)]
  print(paste0("Writing filtered vcf to file ",VCF_OUTPUT_FILE.filtered))
  write.table(VCF.filtered, file=VCF_OUTPUT_FILE.filtered, quote = F, row.names = F, col.names = T, sep = "\t")  
} else {
  if (! SKIP_PLOTS) {
    print("Plotting unfiltered data")
    print("BQD")
    plot_BQD_to_pdf(PDF_OUTPUT_FILE = paste0(PDF_OUTPUT_FILE_PREFIX,".BQD_combined.pdf"), data.bq.triplet = data.bq.triplet, whatToPlot = "BQD")
    print("BQD_CoV")
    plot_BQD_to_pdf(PDF_OUTPUT_FILE = paste0(PDF_OUTPUT_FILE_PREFIX,".BQD_CoV.pdf"), data.bq.triplet = data.bq.triplet, whatToPlot = "BQ_CoV")
    print("BQD_sampleIndividual")
    plot_BQD_to_pdf(PDF_OUTPUT_FILE = paste0(PDF_OUTPUT_FILE_PREFIX,".BQD_individual.pdf"), data.bq.triplet = data.bq.triplet, whatToPlot = "BQD_sampleIndividual")
    print("BQD_sampleIndividual_ColoredByChromosome")
    plot_BQD_to_pdf(PDF_OUTPUT_FILE = paste0(PDF_OUTPUT_FILE_PREFIX,".BQD_individual_CHROMcolored.pdf"), data.bq.triplet = data.bq.triplet, whatToPlot = "BQD_sampleIndividual_ColoredByChromosome")
    print("BQD_sampleIndividual_ColoredByVAF")
    plot_BQD_to_pdf(PDF_OUTPUT_FILE = paste0(PDF_OUTPUT_FILE_PREFIX,".BQD_individual_VAFcolored.pdf"), data.bq.triplet = data.bq.triplet, whatToPlot = "BQD_sampleIndividual_ColoredByVAF")
    print("BQD_sampleIndividual_ColoredByReadPosition")
    ReadPositionQuantile=0.5
    plot_BQD_to_pdf(PDF_OUTPUT_FILE = paste0(PDF_OUTPUT_FILE_PREFIX,".BQD_individual_ReadPosColored_Q",ReadPositionQuantile*100,".pdf"), data.bq.triplet = data.bq.triplet,
                    whatToPlot = "BQD_sampleIndividual_ColoredByReadPosition", ReadPositionQuantile=ReadPositionQuantile )

    ReadPositionQuantile=0.6
    plot_BQD_to_pdf(PDF_OUTPUT_FILE = paste0(PDF_OUTPUT_FILE_PREFIX,".BQD_individual_ReadPosColored_Q",ReadPositionQuantile*100,".pdf"), data.bq.triplet = data.bq.triplet,
                    whatToPlot = "BQD_sampleIndividual_ColoredByReadPosition", ReadPositionQuantile=ReadPositionQuantile )
    
    ReadPositionQuantile=0.7
    plot_BQD_to_pdf(PDF_OUTPUT_FILE = paste0(PDF_OUTPUT_FILE_PREFIX,".BQD_individual_ReadPosColored_Q",ReadPositionQuantile*100,".pdf"), data.bq.triplet = data.bq.triplet,
                    whatToPlot = "BQD_sampleIndividual_ColoredByReadPosition", ReadPositionQuantile=ReadPositionQuantile )
    
    ReadPositionQuantile=0.8
    plot_BQD_to_pdf(PDF_OUTPUT_FILE = paste0(PDF_OUTPUT_FILE_PREFIX,".BQD_individual_ReadPosColored_Q",ReadPositionQuantile*100,".pdf"), data.bq.triplet = data.bq.triplet,
                    whatToPlot = "BQD_sampleIndividual_ColoredByReadPosition", ReadPositionQuantile=ReadPositionQuantile )
    
    print("BQ vs. ReadPosition Scatterplot")
    plot_altBQ_vs_altReadPos(data.bq.triplet = data.bq.triplet, 
                             transitions = transitions, AUC_threshold=0.25, 
                             RESULT_PDF = paste0(MPILEUP_FOLDER,"Scatterplot_BQ_ReadPos.pdf"))
    
  }
}
print("Finished plotting")
