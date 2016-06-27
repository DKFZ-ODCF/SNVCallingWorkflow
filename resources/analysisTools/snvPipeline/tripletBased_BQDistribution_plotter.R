library(getopt)
library(ggplot2)
library(Biostrings) # for reverseComplement


MAX_BASE_QUALITY=55
ALT.MEDIAN.THRESHOLD = 20
parameter.density.bandwidth = 2.00


opt = getopt(matrix(c(
  'vcfInputFile', 'v', 1, "character",
  # 'baseFolder', 'f', 1, "character",
  'mpileupFolder', 'm', 1, "character",
  'alignmentFolder', 'a', 1, "character",
  'PID', 'p', 1, "character",
  'background', 'b', 2, "integer",
  'outFile', 'o', 1, "character",
  'forceRerun', 'R', 2, "integer",
  'combineRevcomp', 'c', 2, "integer"
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
if (is.null(opt$outFile)){      # no vcf file specified
  cat("Please specify the output pdf file.\n"); 
  q(status=2);      # quit, status unequal 0 means error
}


if (COMBINE_REVCOMP) {
  transitions = data.frame( c(rep("C",3),rep("T",3)), c("A","G","T","A","C","G"), stringsAsFactors = F)  
} else {
  transitions = data.frame( c(rep("A",3),rep("C",3),rep("G",3),rep("T",3)), c("C","G","T","A","G","T","A","C","T","A","C","G"), stringsAsFactors = F)
}
colnames(transitions) = c("FROM", "TO")



# RPP_FOLDER="/icgc/dkfzlsdf/analysis/B080/warsow/Cavathrombus/results_per_pid/BW/"
# opt$PID = "BW"
# opt$mpileupFolder = paste0(RPP_FOLDER,"SNV_Calling.backup/mpileup2.3_VRenalisAdjacent_vs_NormalKidney/")
# opt$vcfInputFile = paste0(opt$mpileupFolder,"/snvs_",opt$PID,"_somatic_snvs_conf_8_to_10.vcf")
# opt$alignmentFolder = paste0(RPP_FOLDER,"/alignment/")
# opt$outFile = paste0(opt$mpileupFolder,"snvs_",opt$PID,"_triplet-specific_BQ_distributions.pdf")
# opt$background = F
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
DATA_RESULTS_FILE=paste0(MPILEUP_FOLDER,"BQ_TripletSpecific.RData")

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
  
  #--------------------------------------------------------------------------------------------------------CHR----POS----REF---ALT---------------------BaseBefore-----BaseAfter
  data.bq.triplet = read.table(pipe(paste0("cat ",vcfInputFile, " | cut -f 1,2,4,5,11 | perl -ne '$_ =~ /^(.*?)\t(.*?)\t(.*?)\t(.*?)\t[ACGTNacgtn]{9}([ACGTNacgtn]),([ACGTNacgtn])[ACGTNacgtn]{9}/; print \"$1\t$2\t$3\t$4\t$5\t$6\n\";' ")), comment.char = '', sep = "\t", header = T, stringsAsFactors = F)
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
    
    
    line = cbind(line, countsLine.ref, countsLine.alt)
    
    return (line)
  })
  
  
  data.bq.triplet = do.call(rbind, data.bq.triplet.processed)
  
  remove(data.bq.triplet.processed, data.bq.ref, data.bq.alt)
  
  data.bq.triplet$"BaseBefore" = factor(data.bq.triplet$"BaseBefore", levels=c("A","C","G","T"))
  data.bq.triplet$"REF" = factor(data.bq.triplet$"REF", levels=c("A","C","G","T"))
  data.bq.triplet$"ALT" = factor(data.bq.triplet$"ALT", levels=c("A","C","G","T"))
  data.bq.triplet$"BaseAfter" = factor(data.bq.triplet$"BaseAfter", levels=c("A","C","G","T"))
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
data.bq.triplet.filtered = do.call(rbind, apply(transitions, 1, function(TRANSITION) {
  if ( "THIS SNV IS ARTIFACTED" == "THIS SNV IS ARTIFACTED" ) {
    # TODO: only filter channels that show artifact (criterion still has to be defined)
    SNVs = data.bq.triplet[with(data.bq.triplet, substr(REVCOMP_TRIPLET,2,2)==TRANSITION[1] & substr(REVCOMP_TRIPLET,3,3)==TRANSITION[2]),]
    indices = which(SNVs$MeadianBQ.alt >= ALT.MEDIAN.THRESHOLD)
    SNVs = SNVs[indices,]    
    return (SNVs)
  }
}))
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


plot_BQD_to_pdf = function(PDF_OUTPUT_FILE, data.bq.triplet) {
  pdf(PDF_OUTPUT_FILE, 8.27, 6.0)
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
  dev.off()
}



plot_BQD_to_pdf(PDF_OUTPUT_FILE = PDF_OUTPUT_FILE, data.bq.triplet = data.bq.triplet)

PDF_OUTPUT_FILE.filtered = sub(".pdf", paste0("_filteredAltMedian",ALT.MEDIAN.THRESHOLD,".pdf"), PDF_OUTPUT_FILE)
plot_BQD_to_pdf(PDF_OUTPUT_FILE = PDF_OUTPUT_FILE.filtered, data.bq.triplet = data.bq.triplet.filtered)


VCF_OUTPUT_FILE.filtered = sub(".vcf", paste0("_filteredAltMedian",ALT.MEDIAN.THRESHOLD,".vcf"), vcfInputFile)
VCF = read.table(vcfInputFile, comment.char = '', sep = "\t", header = T, stringsAsFactors = F, check.names = F)
colnames(VCF)[1] = "CHROM"
VCF.filtered = merge(VCF, data.bq.triplet.filtered[,c("CHROM","POS","REF","ALT")])
chrOrder <-c((1:22),"X","Y")
VCF.filtered$CHROM =  factor(VCF.filtered$CHROM, chrOrder, ordered=TRUE)
VCF.filtered = VCF.filtered[with(data=VCF.filtered, order(CHROM, POS)),]
colnames(VCF.filtered)[1] = "#CHROM"
ncol=ncol(VCF.filtered)
VCF.filtered = VCF.filtered[,c(1,2,5,3,4,6:ncol)]
write.table(VCF.filtered, file=VCF_OUTPUT_FILE.filtered, quote = F, row.names = F, col.names = T, sep = "\t")

