library(getopt)
# library(ggplot2)
library(Biostrings) # for reverseComplement
# library(grid) # for unit,gpar
library(reshape2)

options(stringsAsFactors = FALSE)
options(scipen=999)


SEQUENCE_CONTEXT_COLUMN_INDEX=11
MAX_BASE_QUALITY=55


opt = getopt(matrix(c(
  'vcfInputFile', 'v', 1, "character",
  'mpileupFolder', 'm', 1, "character",
  'PID', 'p', 1, "character",
  'outFile', 'o', 1, "character",
  'combineRevcomp', 'c', 2, "integer",
  'sequenceContextColumnIndex', 's', 2, "integer"
),ncol=4,byrow=TRUE));


# RPP_FOLDER="/icgc/pcawg/analysis/train2full/projects/CMDI-UK/"
# opt$PID = "f8e61a02-654c-c226-e040-11ac0d481b60"
# opt$PID = "f9012277-55c4-2486-e040-11ac0c48236f"
# opt$mpileupFolder = paste0(RPP_FOLDER,opt$PID,"/merged_may2016Calls_filtered_BQD_plots")
# opt$vcfInputFile = paste0(opt$mpileupFolder,"/",opt$PID,"_somatic.snv_mnv.vcf_filtered.OnlyPassed.vcf")
# opt$outFile = paste0(opt$mpileupFolder,"/snvs_",opt$PID,"_SangerArtifcatDetected.txt")
# opt$combineRevcomp = 1

# RPP_FOLDER="/icgc/pcawg/analysis/train2full/projects/BLCA-US/"
# opt$PID = "2b142863-b963-4cc9-8f8f-c72503c93390"
# opt$PID = "0c7aca3f-e006-4de3-afc2-20b4f727d4fd"
# opt$mpileupFolder = paste0(RPP_FOLDER,opt$PID,"/merged_may2016Calls_filtered_BQD_plots")
# opt$vcfInputFile = paste0(opt$mpileupFolder,"/",opt$PID,"_somatic.snv_mnv.vcf_filtered.OnlyPassed.vcf")
# opt$outFile = paste0(opt$mpileupFolder,"/snvs_",opt$PID,"_SangerArtifcatDetected.txt")
# opt$combineRevcomp = 1
# opt$sequenceContextColumnIndex = 10

# RPP_FOLDER="/icgc/pcawg/analysis/train2full/projects/CLLE-ES/"
# opt$PID = "873ed251-fa44-465f-99e0-c90a42ea13d9"
# opt$mpileupFolder = paste0(RPP_FOLDER,opt$PID,"/merged_may2016Calls_filtered_BQD_plots")
# opt$vcfInputFile = paste0(opt$mpileupFolder,"/",opt$PID,"_somatic.snv_mnv.vcf_filtered.OnlyPassed.vcf")
# opt$outFile = paste0(opt$mpileupFolder,"/snvs_",opt$PID,"_SangerArtifcatDetected.txt")
# opt$combineRevcomp = 1

# RPP_FOLDER="/icgc/pcawg/analysis/train2full/projects/BOCA-UK/"
# opt$PID = "f85ae2b7-cebf-17a2-e040-11ac0c48033a"
# opt$PID = "fc970467-ae01-4a8f-e040-11ac0c48163b"
# opt$PID = "f82d213f-bc99-5b1d-e040-11ac0c486880"
# opt$mpileupFolder = paste0(RPP_FOLDER,opt$PID,"/merged_may2016Calls_filtered_BQD_plots")
# opt$vcfInputFile = paste0(opt$mpileupFolder,"/",opt$PID,"_somatic.snv_mnv.vcf_filtered.OnlyPassed.vcf")
# opt$vcfInputFile = paste0(opt$mpileupFolder,"/",opt$PID,"_somatic.snv_mnv.vcf_filtered.OnlyPassed.withoutHeader.vcf")
# opt$outFile = paste0(opt$mpileupFolder,"/snvs_",opt$PID,"_SangerArtifcatDetected.txt")
# opt$combineRevcomp = 1

# RPP_FOLDER="/icgc/pcawg/analysis/train2full/projects/LAML-US/"
# opt$PID = "fc2d6c31-dadb-45cb-9053-c87f535a77a5"
# opt$PID = "fc970467-ae01-4a8f-e040-11ac0c48163b"
# opt$PID = "f82d213f-bc99-5b1d-e040-11ac0c486880"
# opt$mpileupFolder = "/icgc/pcawg/analysis/train2full/projects/LAML-US/fc2d6c31-dadb-45cb-9053-c87f535a77a5/artifactDetection/Sanger"
# opt$vcfInputFile = "/icgc/pcawg/analysis/train2full/projects/LAML-US/fc2d6c31-dadb-45cb-9053-c87f535a77a5/artifactDetection/Sanger/snvs_fc2d6c31-dadb-45cb-9053-c87f535a77a5_somatic_snvs_conf_8_to_10.OnlyPassed.withoutHeader.vcf"
# opt$outFile = "/icgc/pcawg/analysis/train2full/projects/LAML-US/fc2d6c31-dadb-45cb-9053-c87f535a77a5/artifactDetection/Sanger/snvs_fc2d6c31-dadb-45cb-9053-c87f535a77a5_SangerArtifcatData.txt"
# opt$combineRevcomp = 1
# opt$sequenceContextColumnIndex = 10

if (is.null(opt$vcfInputFile)){
  cat("Please specify the file that contains the SNVs for which the base score distribution plot shall be created.\n"); 
  q(status=1);      # quit, status unequal 0 means error
}
if (is.null(opt$mpileupFolder)){
  cat("Please specify the mpileup folder.\n"); 
  q(status=1);      # quit, status unequal 0 means error
}
if (is.null(opt$PID)){
  cat("Please specify the PID.\n"); 
  q(status=1);      # quit, status unequal 0 means error
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
if (is.null(opt$outFile)){      # no vcf file specified
  cat("Please specify the output file.\n"); 
  q(status=2);      # quit, status unequal 0 means error
}
if (! is.null(opt$sequenceContextColumnIndex)){
  tmp = as.integer(opt$sequenceContextColumnIndex)
  if (!is.na(tmp)) {
    SEQUENCE_CONTEXT_COLUMN_INDEX = tmp
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



PID=opt$PID
MPILEUP_FOLDER=paste0(opt$mpileupFolder,"/")
vcfInputFile = paste0(opt$vcfInputFile)
RefAlleleBaseQualitiesFile=paste0(MPILEUP_FOLDER,"snvs_",PID,"_reference_allele_base_qualities.txt.gz", collapse = "")
AltAlleleBaseQualitiesFile=paste0(MPILEUP_FOLDER,"snvs_",PID,"_alternative_allele_base_qualities.txt.gz", collapse = "")
OUTPUT_FILE=paste0(opt$outFile)



DATA_RESULTS_FILE=paste0(MPILEUP_FOLDER,"BQ_TripletSpecific.RData")
if ( file.exists(DATA_RESULTS_FILE) ) {
  lnames = load(file = DATA_RESULTS_FILE)
} else {
  
  # READ IN DATA FOR CHANNEL-SPECIFIC BQs
  #--------------------------------------------------------------------------------------------------------CHR----POS----REF---ALT---------------------BaseBefore-----BaseAfter
  data.bq.triplet = read.table(pipe(paste0("cat ",vcfInputFile, " | grep -v '^##' | cut -f 1,2,4,5,",SEQUENCE_CONTEXT_COLUMN_INDEX," | perl -ne '$_ =~ /^(.*?)\t(.*?)\t(.*?)\t(.*?)\t[ACGTNacgtn]{9}([ACGTNacgtn]),([ACGTNacgtn])[ACGTNacgtn]{9}/; print \"$1\t$2\t$3\t$4\t$5\t$6\n\";' ")), comment.char = '', sep = "\t", header = T, stringsAsFactors = F)
  if (all(is.na(data.bq.triplet))) {
    cat("Please check the parameter -s [sequenceContextColumnIndex]!")
    quit(save = F, status = 11)
  }
  colnames(data.bq.triplet) = c("CHROM", "POS","REF","ALT","BaseBefore","BaseAfter")
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
  
  remove(data.bq.ref, data.bq.alt)
  
  
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
  
  remove(data.bq.triplet.processed)
  
  data.bq.triplet$"BaseBefore" = factor(data.bq.triplet$"BaseBefore", levels=c("A","C","G","T"))
  data.bq.triplet$"REF" = factor(data.bq.triplet$"REF", levels=c("A","C","G","T"))
  data.bq.triplet$"ALT" = factor(data.bq.triplet$"ALT", levels=c("A","C","G","T"))
  data.bq.triplet$"BaseAfter" = factor(data.bq.triplet$"BaseAfter", levels=c("A","C","G","T"))
  data.bq.triplet$"Triplet" = as.factor(as.character(data.bq.triplet$"Triplet"))
  
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


AUC.alt.DF = data.frame(Triplet = character(0), 
                             AUC30=numeric(0), AUC5.30=numeric(0), AUC13.30=numeric(0), 
                             AUC25=numeric(0), AUC5.25=numeric(0), AUC13.25=numeric(0), 
                             AUC20=numeric(0), AUC5.20=numeric(0), AUC13.20=numeric(0) )

data.SangerArtifact = apply(transitions, 1, function(transition) {
  # transition=transitions[i,]
  # transition=transitions[6,]
  from=transition["FROM"]
  to=transition["TO"]
  fromto = paste0(from,to)
  n.total = sum(substr(data.bq.triplet$REVCOMP_TRIPLET,2,2)==as.character(from) & substr(data.bq.triplet$REVCOMP_TRIPLET,3,3)==as.character(to))
  data.SangerArtifact.tmp = data.frame(Triplet=character(0), count=integer(0), proportion=numeric(0), AUC.BQ30.ref=numeric(0), AUC.BQ30.alt=numeric(0))
  for (baseBefore in c("A","C","G","T")) {
    for (baseAfter in c("A","C","G","T")) {
      # baseBefore="G"
      # baseAfter="G"
      triplet = paste0(baseBefore,from,to,baseAfter, collapse = "")
      transitionSubset = data.bq.triplet[data.bq.triplet$REVCOMP_TRIPLET == triplet,]
      n = nrow(transitionSubset)

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

      proportion.channel = round(n/n.total,3)
      
      AUC.BQ30.ref = round(sum(molten[with(molten, BQ <= 30 & type == 'REF'),"value"]),3)
      AUC.BQ30.alt = round(sum(molten[with(molten, BQ <= 30 & type == 'ALT'),"value"]),3)
      AUC.BQ5.30.ref = round(sum(molten[with(molten, BQ >= 5 & BQ <= 30 & type == 'REF'),"value"]),3)
      AUC.BQ5.30.alt = round(sum(molten[with(molten, BQ >= 5 & BQ <= 30 & type == 'ALT'),"value"]),3)
      AUC.BQ13.30.ref = round(sum(molten[with(molten, BQ >= 13 & BQ <= 30 & type == 'REF'),"value"]),3)
      AUC.BQ13.30.alt = round(sum(molten[with(molten, BQ >= 13 & BQ <= 30 & type == 'ALT'),"value"]),3)            
      
      AUC.BQ25.ref = round(sum(molten[with(molten, BQ <= 25 & type == 'REF'),"value"]),3)
      AUC.BQ25.alt = round(sum(molten[with(molten, BQ <= 25 & type == 'ALT'),"value"]),3)
      AUC.BQ5.25.ref = round(sum(molten[with(molten, BQ >= 5 & BQ <= 25 & type == 'REF'),"value"]),3)
      AUC.BQ5.25.alt = round(sum(molten[with(molten, BQ >= 5 & BQ <= 25 & type == 'ALT'),"value"]),3)
      AUC.BQ13.25.ref = round(sum(molten[with(molten, BQ >= 13 & BQ <= 25 & type == 'REF'),"value"]),3)
      AUC.BQ13.25.alt = round(sum(molten[with(molten, BQ >= 13 & BQ <= 25 & type == 'ALT'),"value"]),3)                  
      
      AUC.BQ20.ref = round(sum(molten[with(molten, BQ <= 20 & type == 'REF'),"value"]),3)
      AUC.BQ20.alt = round(sum(molten[with(molten, BQ <= 20 & type == 'ALT'),"value"]),3)
      AUC.BQ5.20.ref = round(sum(molten[with(molten, BQ >= 5 & BQ <= 20 & type == 'REF'),"value"]),3)
      AUC.BQ5.20.alt = round(sum(molten[with(molten, BQ >= 5 & BQ <= 20 & type == 'ALT'),"value"]),3)      
      AUC.BQ13.20.ref = round(sum(molten[with(molten, BQ >= 13 & BQ <= 20 & type == 'REF'),"value"]),3)
      AUC.BQ13.20.alt = round(sum(molten[with(molten, BQ >= 13 & BQ <= 20 & type == 'ALT'),"value"]),3)                        
      
      
      AUC.alt.DF <<- rbind(AUC.alt.DF, c(triplet, 
                                                   AUC.BQ30.alt, AUC.BQ5.30.alt, AUC.BQ13.30.alt, 
                                                   AUC.BQ25.alt, AUC.BQ5.25.alt, AUC.BQ13.25.alt, 
                                                   AUC.BQ20.alt, AUC.BQ5.20.alt, AUC.BQ13.20.alt))

      data.SangerArtifact.tmp = rbind(data.SangerArtifact.tmp, t(c(triplet, n, proportion.channel, 
                                                                   AUC.BQ30.ref, AUC.BQ30.alt, AUC.BQ5.30.ref, AUC.BQ5.30.alt, AUC.BQ13.30.ref, AUC.BQ13.30.alt,
                                                                   AUC.BQ25.ref, AUC.BQ25.alt, AUC.BQ5.25.ref, AUC.BQ5.25.alt, AUC.BQ13.25.ref, AUC.BQ13.25.alt,
                                                                   AUC.BQ20.ref, AUC.BQ20.alt, AUC.BQ5.20.ref, AUC.BQ5.20.alt, AUC.BQ13.20.ref, AUC.BQ13.20.alt)))
    }
  }
  return (data.SangerArtifact.tmp)
})


data.SangerArtifact = do.call(rbind, data.SangerArtifact)
colnames(data.SangerArtifact) = c("Triplet", "count", "proportion", 
                                  "AUC.BQ30.ref", "AUC.BQ30.alt", "AUC.BQ5.30.ref", "AUC.BQ5.30.alt", "AUC.BQ13.30.ref", "AUC.BQ13.30.alt", 
                                  "AUC.BQ25.ref", "AUC.BQ25.alt", "AUC.BQ5.25.ref", "AUC.BQ5.25.alt", "AUC.BQ13.25.ref", "AUC.BQ13.25.alt", 
                                  "AUC.BQ20.ref", "AUC.BQ20.alt", "AUC.BQ5.20.ref", "AUC.BQ5.20.alt", "AUC.BQ13.20.ref", "AUC.BQ13.20.alt")
rownames(data.SangerArtifact) = data.SangerArtifact$Triplet

data.SangerArtifact$count = as.integer(data.SangerArtifact$count)
for (colname in colnames(data.SangerArtifact)[c(-1,-2)]) {
  data.SangerArtifact[,colname] = ifelse(data.SangerArtifact[,colname]=="0", NA, as.numeric(data.SangerArtifact[,colname]))
}


colnames(AUC.alt.DF) = c("Triplet", "AUC30", "AUC5.30", "AUC13.30", "AUC25", "AUC5.25", "AUC13.25", "AUC20", "AUC5.20", "AUC13.20")
rownames(AUC.alt.DF) = AUC.alt.DF$Triplet
for (colname in colnames(AUC.alt.DF)[c(-1)]) {
  AUC.alt.DF[,colname] = as.numeric(AUC.alt.DF[,colname])
  AUC.alt.DF[,colname] = ifelse(AUC.alt.DF[,colname]=="0", NA, as.numeric(AUC.alt.DF[,colname]))
}
AUC.alt.DF$const = 1

pdf(file=paste0(MPILEUP_FOLDER,"snvs_",PID,"_zScoreBased_AUC.BQ30_BiasPlot.pdf", collapse = ""), paper = "a4r")
# for (upperThresholdBQ in c(20,25,30)) {
#   for (lowerThresholdBQ in c(0,5,13)) {
for (upperThresholdBQ in c(30)) {
  for (lowerThresholdBQ in c(0)) {    
    if (lowerThresholdBQ > 0 ) {
      # lmCommand = paste0("mod = lm(AUC",lowerThresholdBQ,".",upperThresholdBQ," ~ const, data=AUC.alt.DF)")
      # CooksD_colName = paste0("CooksD",lowerThresholdBQ,".",upperThresholdBQ)
      # FoldCookMean_colName = paste0("FoldCookMean",lowerThresholdBQ,".",upperThresholdBQ)      
      # AUC_colName = paste0("AUC.BQ",lowerThresholdBQ,".",upperThresholdBQ,".alt")
      # mainTitleForPlot = paste0("Cook's Distance for AUC.BQ",lowerThresholdBQ,".",upperThresholdBQ,".alt\nfor PID ",PID)
    } else {
      # lmCommand = paste0("mod = lm(AUC",upperThresholdBQ," ~ const, data=AUC.alt.DF)")
      # CooksD_colName = paste0("CooksD",upperThresholdBQ)
      FoldCookMean_colName = paste0("FoldCookMean",upperThresholdBQ)
      Residuum_colName = paste0("Residuum.BQ",upperThresholdBQ)
      FoldResiduumMean_colname = paste0("FoldResiduumMean.BQ",upperThresholdBQ)
      zScore_colname = paste0("zScore.BQ",upperThresholdBQ)
      AUC_colName = paste0("AUC.BQ",upperThresholdBQ,".alt")
      # mainTitleForPlot = paste0("Cook's Distance for AUC.BQ",upperThresholdBQ,".alt\nfor PID ",PID)
      mainTitleForPlot = paste0("AUC.BQ",upperThresholdBQ,".alt with zScore-based outlier detection\nfor PID ",PID)
    }
    # eval(parse(text=lmCommand))
    # cooksd <- cooks.distance(mod)
    # http://r-statistics.co/Outlier-Treatment-With-R.html
    AUC = data.SangerArtifact[,AUC_colName]
    AUC.mean = mean(AUC, na.rm=T)
    AUC.stdv = sd(AUC,na.rm = T)
    AUC.zScore = abs(AUC - AUC.mean) / AUC.stdv
    AUC.alt.DF[,zScore_colname] = AUC.zScore
    # AUC.alt.DF[names(cooksd),CooksD_colName] = cooksd
    # data.SangerArtifact[,FoldCookMean_colName] = round(AUC.alt.DF[,CooksD_colName]/mean(AUC.alt.DF[,CooksD_colName], na.rm=T),2)
    
    # indices_orange = cooksd>4*mean(cooksd, na.rm=T)
    # indices_orange = names(indices_orange)[indices_orange]
    # indices_red = cooksd>10*mean(cooksd, na.rm=T)
    # indices_red = names(indices_red)[indices_red]
    
    indices_orange = AUC.alt.DF[,zScore_colname]>2
    indices_orange = rownames(AUC.alt.DF)[indices_orange]
    indices_orange = indices_orange[!is.na(indices_orange)]
    indices_red = AUC.alt.DF[,zScore_colname]>4
    indices_red = rownames(AUC.alt.DF)[indices_red]    
    indices_red = indices_red[!is.na(indices_red)]
    
    AUC.alt.DF$printLabel = ""
    AUC.alt.DF$printAUC = ""
    AUC.alt.DF$color = "black"
    if (length(indices_orange)>0) {
      AUC.alt.DF[indices_orange, "printLabel"] = AUC.alt.DF[indices_orange,"Triplet"]
      AUC.alt.DF[indices_orange, "printAUC"] = paste0("AUC: ",data.SangerArtifact[indices_orange,AUC_colName],"  n: ",data.SangerArtifact[indices_orange,"count"])
      AUC.alt.DF[indices_orange, "color"] = "orange"      
    }
    if (length(indices_red)) {
      AUC.alt.DF[indices_red, "color"] = "red"
    }
    
        
    # plot(AUC.alt.DF[,CooksD_colName], pch="*", cex=2, col=AUC.alt.DF$color, ylab="CooksDistance", xlab="Triplet Index", ylim=c(0.0,max(c(0.6,max(cooksd)+0.1))), main=mainTitleForPlot)
    # abline(h = 4*mean(cooksd, na.rm=T), col="orange")
    # text(96, 4*mean(cooksd, na.rm=T), "4-fold", col="orange", cex=0.6)
    # abline(h = 10*mean(cooksd, na.rm=T), col="red")
    # text(95, 10*mean(cooksd, na.rm=T), "10-fold", col="red", cex=0.6)
    # for (fold in c(25,50,75,100)) {
    #   abline(h = fold*mean(cooksd, na.rm=T), col="grey")
    #   text(95, fold*mean(cooksd, na.rm=T), paste0(fold,"-fold"), col="grey", cex=0.6)      
    # }
    # abline(h = mean(cooksd, na.rm=T), col="grey")
    # text(95, mean(cooksd, na.rm=T), paste0("mean"), col="grey", cex=0.6)      
    # text(x=1:nrow(AUC.alt.DF), y=AUC.alt.DF[,CooksD_colName]+0.030, labels=AUC.alt.DF$printLabel, col=AUC.alt.DF$color, cex=0.8)
    # text(x=1:nrow(AUC.alt.DF), y=AUC.alt.DF[,CooksD_colName]+0.015, labels=AUC.alt.DF$printAUC, col=AUC.alt.DF$color, cex=0.6)
    
    
    plot(data.SangerArtifact[,AUC_colName], pch="*", cex=2, col=AUC.alt.DF$color, 
         ylab="AUC", xlab="Triplet Index", ylim=c(0,1.05), 
         main=mainTitleForPlot)    
    abline(h = 2*AUC.stdv+AUC.mean, col="orange")
    text(96, 2*AUC.stdv+AUC.mean+0.01, "zScore 2", col="orange", cex=0.6)
    abline(h = 4*AUC.stdv+AUC.mean, col="red")
    text(96, 4*AUC.stdv+AUC.mean+0.01, "zScore 4", col="red", cex=0.6)

    abline(h = AUC.mean, col="grey")
    text(95, AUC.mean+0.01, paste0("mean AUC"), col="grey", cex=0.6)      
    text(x=1:nrow(AUC.alt.DF), y=data.SangerArtifact[,AUC_colName]+0.04, labels=AUC.alt.DF$printLabel, col=AUC.alt.DF$color, cex=0.8)
    text(x=1:nrow(AUC.alt.DF), y=data.SangerArtifact[,AUC_colName]+0.021, labels=AUC.alt.DF$printAUC, col=AUC.alt.DF$color, cex=0.6)    
  }
}
dev.off()




data.SangerArtifact = merge(data.SangerArtifact, AUC.alt.DF, by="Triplet")

data.SangerArtifact = data.SangerArtifact[order(substr(data.SangerArtifact$Triplet,start = 2,3),substr(data.SangerArtifact$Triplet,start = 1,1),substr(data.SangerArtifact$Triplet,start = 4,4)),]
write.table(data.SangerArtifact, file = OUTPUT_FILE, sep = "\t", row.names = F, col.names = T, quote = F)



# if (data.SangerArtifact[data.SangerArtifact$Triplet=="GTGG","proportion"] > 0.2) {
#   if (data.SangerArtifact[data.SangerArtifact$Triplet=="GTGG","count"] > 5) {
#     if (data.SangerArtifact[data.SangerArtifact$Triplet=="GTGG","AUC.BQ30.alt"] > 0.40) {
#       data.SangerArtifact = data.SangerArtifact[order(substr(data.SangerArtifact$Triplet,start = 2,3),substr(data.SangerArtifact$Triplet,start = 1,1),substr(data.SangerArtifact$Triplet,start = 4,4)),]
#       write.table(data.SangerArtifact, file = OUTPUT_FILE, sep = "\t", row.names = F, col.names = T, quote = F)
#     }
#   }
# }

# if (data.SangerArtifact[data.SangerArtifact$Triplet=="GTGG","proportion"] > 0.0) {
#   if (data.SangerArtifact[data.SangerArtifact$Triplet=="GTGG","count"] > 0) {
#     if (data.SangerArtifact[data.SangerArtifact$Triplet=="GTGG","AUC.BQ30.alt"] > 0.00) {
#       data.SangerArtifact = data.SangerArtifact[order(substr(data.SangerArtifact$Triplet,start = 2,3),substr(data.SangerArtifact$Triplet,start = 1,1),substr(data.SangerArtifact$Triplet,start = 4,4)),]
#       write.table(data.SangerArtifact, file = OUTPUT_FILE, sep = "\t", row.names = F, col.names = T, quote = F)
#     }
#   }          
# }




