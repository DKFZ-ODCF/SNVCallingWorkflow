#!/usr/bin/env Rscript
#
# Copyright (c) 2018 German Cancer Research Center (DKFZ).
#
# Distributed under the GPL-2/GPL-3 License (http://www.gnu.de/documents/gpl-2.0.en.html, https://www.gnu.org/licenses/gpl-3.0.en.html).
#

library(getopt)
library(ggplot2)
library(gridExtra) # for tableGrob
library(grid) # for unit,gpar
library(data.table) # for rbindlist

library(Biostrings) # for reverseComplement

combineReverseComplement = T
Plottype.Ratio.UniColor=1
Plottype.Differences.BiColor=2
SEQUENCE_CONTEXT_COLUMN_INDEX=11


## get parameters
opt = getopt(matrix(c(
  # 'pid', 'p', 1, "character",
  'vcfInputFile', 'v', 1, "character",
  'refScores', 'r', 1, "character",
  'altScores', 'a', 1, "character",
  'baseScoreThreshold', 't', 2, "integer",
  'plotType', 'p', 2, "character",
  'descriptionForMainTitle', 'd', 2, "character",
  'outFile', 'o', 1, "character",
  'sequenceContextColumnIndex', 's', 2, "integer"
),ncol=4,byrow=TRUE));



checkForMissingParameter = function(parameter, errorText, exitCode=1) {
  if (is.null(opt[[parameter]])){
    cat(paste0(errorText,"\n")) 
    q(exitCode)      # quit, status unequal 0 means error
  }
}

setIntegerValueFromParameter = function(parameterName, valueName) {
  if (! is.null(opt[[parameterName]])){
    tmp = as.integer(opt[[parameterName]])
    if (!is.na(tmp)) {
      assign(valueName, tmp)
    }
  }  
}

checkForMissingParameter("vcfInputFile", "Please specify the file that contains the SNVs for which the base score distribution plot shall be created.", 1)
checkForMissingParameter("refScores", "Please specify the file containing reference allele base scores.", 1)
checkForMissingParameter("altScores", "Please specify the file containing alternative allele base scores.", 1)
checkForMissingParameter("outFile", "Please specify the output pdf file.", 1)
checkForMissingParameter("baseScoreThreshold", "No base score threshold is set", 1)


if (is.null(opt$plotType)){      # noplot type is set
  plotType=Plottype.Differences.BiColor
} else if (opt$plotType == "Differences") {
  plotType=Plottype.Differences.BiColor
} else if (opt$plotType == "Ratios") {
  plotType=Plottype.Ratio.UniColor
} else {
  cat("Please specify the plotType as one of \"Differences\" or \"Ratios\".\n"); 
  q(status=3);      # quit, status unequal 0 means error
}

setIntegerValueFromParameter(parameterName = "sequenceContextColumnIndex", valueName = "SEQUENCE_CONTEXT_COLUMN_INDEX")


options(stringsAsFactors = FALSE)


data = try( read.table(pipe(paste0("cut -f 1,2,4,5,",SEQUENCE_CONTEXT_COLUMN_INDEX," ",opt$vcfInputFile)), comment.char = '', sep = "\t", header = T) )
colnames(data)[1] = "CHROM"


refScores = read.table(opt$refScores)
colnames(refScores) = c("CHROM", "POS", "BaseScoreString.ref")
altScores = read.table(opt$altScores)
colnames(altScores) = c("CHROM", "POS", "BaseScoreString.alt")

# MERGE commands ensure that we only keep the necessary positions (in case we extracted the bases for all SNVs and not only the HC somatic SNVs)
refScores = merge(refScores, data[,1:2], by=c("CHROM","POS"))
altScores = merge(altScores, data[,1:2], by=c("CHROM","POS"))


refScores$BaseScore.ref = apply(refScores, 1, function(line) {
  baseScores = as.integer(unlist(strsplit(unlist(strsplit(line["BaseScoreString.ref"], ";"))[1], ",")))
  baseScores[which(baseScores >= opt$baseScoreThreshold)]
})

refScores$meanBQ.ref = apply(refScores, 1, function(line) { mean(unlist(line["BaseScore.ref"])) })
refScores$medianBQ.ref = apply(refScores, 1, function(line) { median(unlist(line["BaseScore.ref"])) })

altScores$BaseScore.alt = apply(altScores, 1, function(line) {
  baseScores = as.integer(unlist(strsplit(unlist(strsplit(line["BaseScoreString.alt"], ";"))[1], ",")))
  baseScores[which(baseScores >= opt$baseScoreThreshold)]
})
altScores$meanBQ.alt = apply(altScores, 1, function(line) { mean(unlist(line["BaseScore.alt"])) })
altScores$medianBQ.alt = apply(altScores, 1, function(line) { median(unlist(line["BaseScore.alt"])) })


data$baseBefore = apply(data, 1, function(line) {substr(unlist(strsplit(line["SEQUENCE_CONTEXT"], ","))[1], 10, 10)})
data$baseAfter = apply(data, 1, function(line) {substr(unlist(strsplit(line["SEQUENCE_CONTEXT"], ","))[2], 1, 1)})
data$triplet = as.factor(apply(data, 1, function(line) { paste0(line["baseBefore"], line["REF"], line["ALT"], line["baseAfter"], collapse = "") }))

data = merge(data, refScores[,c("CHROM","POS","BaseScoreString.ref","meanBQ.ref","medianBQ.ref")], by=c("CHROM","POS"))
data = merge(data, altScores[,c("CHROM","POS","BaseScoreString.alt","meanBQ.alt","medianBQ.alt")], by=c("CHROM","POS"))

tripletSpecificBaseQualities.ref = aggregate(BaseScoreString.ref ~ triplet, data, c)
rownames(tripletSpecificBaseQualities.ref) = tripletSpecificBaseQualities.ref$triplet
tripletSpecificBaseQualities.alt = aggregate(BaseScoreString.alt ~ triplet, data, c)
rownames(tripletSpecificBaseQualities.alt) = tripletSpecificBaseQualities.alt$triplet

if (any(rownames(tripletSpecificBaseQualities.ref) != rownames(tripletSpecificBaseQualities.alt))) {
  stop("Different rownames for REF an ALT!") 
}


# Combine triplets to reverse complements if requested to do so via the parameters
if (combineReverseComplement) {
  possible_mutations = c("CA", "CG", "CT", "TA", "TC", "TG")
  forbidden_mutations = sapply(possible_mutations, function(mutation) {as.character(complement(DNAString(mutation)))})
  
  allowedTripletsIndices = which(substr(tripletSpecificBaseQualities.ref$triplet, 2,3) %in% possible_mutations)
  tripletSpecificBaseQualityRatios = as.data.frame(tripletSpecificBaseQualities.ref[allowedTripletsIndices,"triplet"])
    colnames(tripletSpecificBaseQualityRatios) = "triplet"
    rownames(tripletSpecificBaseQualityRatios) = tripletSpecificBaseQualityRatios$triplet  
} else {
  possible_mutations = c("AC","AG","AT", "CA","CG","CT", "GA","GC","GT", "TA","TC","TG")
  tripletSpecificBaseQualityRatios = as.data.frame(tripletSpecificBaseQualities.ref$triplet)
    colnames(tripletSpecificBaseQualityRatios) = "triplet"
    rownames(tripletSpecificBaseQualityRatios) = tripletSpecificBaseQualityRatios$triplet  
}




for (triplet in rownames(tripletSpecificBaseQualities.ref)) {
  # do the following for each triplet...
  if (combineReverseComplement) {
    if (substr(triplet,2,3) %in% possible_mutations) {
      # collect the scores for the allowed half of this triplet
      scores.ref = unlist(as.vector(apply(as.data.frame(tripletSpecificBaseQualities.ref[triplet,"BaseScoreString.ref"]), 1, function(BQset) {
        as.integer(unlist(strsplit(unlist(strsplit(BQset, ";"))[1], ",")))
      })))
      scores.alt = unlist(as.vector(apply(as.data.frame(tripletSpecificBaseQualities.alt[triplet,"BaseScoreString.alt"]), 1, function(BQset) {
        as.integer(unlist(strsplit(unlist(strsplit(BQset, ";"))[1], ",")))
      })))
      # collect the scores for the forbidden half of this triplet
      triplet.revComp.tmp = as.character(reverseComplement(DNAString(triplet)))
      triplet.revComp = paste0(substr(triplet.revComp.tmp,1,1),substr(triplet.revComp.tmp,3,3),substr(triplet.revComp.tmp,2,2),substr(triplet.revComp.tmp,4,4))
      tripletSpecificBaseQualities.ref.subset = as.data.frame(tripletSpecificBaseQualities.ref[triplet.revComp,"BaseScoreString.ref"])
      # We have to ensure that the revComp triplet exists. This does not necessarily have to be the case!
      if (nrow(tripletSpecificBaseQualities.ref.subset) > 0) {
        scores.ref = c(scores.ref, unlist(apply(tripletSpecificBaseQualities.ref.subset, 1, function(BQset) {
          as.integer(unlist(strsplit(unlist(strsplit(BQset, ";"))[1], ",")))
        })))  
      }
      tripletSpecificBaseQualities.alt.subset = as.data.frame(tripletSpecificBaseQualities.alt[triplet.revComp,"BaseScoreString.alt"])
      if (nrow(tripletSpecificBaseQualities.alt.subset) > 0) {
        scores.alt = c(scores.alt, unlist(apply(tripletSpecificBaseQualities.alt.subset, 1, function(BQset) {
          as.integer(unlist(strsplit(unlist(strsplit(BQset, ";"))[1], ",")))
        })))        
      }
    } else {
      # this triplet is forbidden and has been handled / will be handled in the corresponding allowed triplet
      next()
    }
  } else {
    scores.ref = unlist(apply(as.data.frame(tripletSpecificBaseQualities.ref[triplet,"BaseScoreString.ref"]), 1, function(BQset) {
      as.integer(unlist(strsplit(unlist(strsplit(BQset, ";"))[1], ",")))
    }))
    scores.alt = unlist(apply(as.data.frame(tripletSpecificBaseQualities.alt[triplet,"BaseScoreString.alt"]), 1, function(BQset) {
      as.integer(unlist(strsplit(unlist(strsplit(BQset, ";"))[1], ",")))
    }))    
  }
  
  # generate some statistics for the current triplet
  mean.ref = mean(scores.ref)
  mean.alt = mean(scores.alt)
  median.ref = median(scores.ref)
  median.alt = median(scores.alt)
  min.ref = min(scores.ref)
  min.alt = min(scores.alt)
  max.ref = max(scores.ref)
  max.alt = max(scores.alt)  
  tripletSpecificBaseQualityRatios[triplet,"min.ref"] = min.ref
  tripletSpecificBaseQualityRatios[triplet,"max.ref"] = max.ref
  tripletSpecificBaseQualityRatios[triplet,"mean.ref"] = mean.ref
  tripletSpecificBaseQualityRatios[triplet,"median.ref"] = median.ref
  tripletSpecificBaseQualityRatios[triplet,"n.ref"] = length(scores.ref)
  tripletSpecificBaseQualityRatios[triplet,"existingBaseScores.ref"] = paste0(names(table(scores.ref)), collapse = ",")
  
  tripletSpecificBaseQualityRatios[triplet,"mean.alt"] = mean.alt
  tripletSpecificBaseQualityRatios[triplet,"min.alt"] = min.alt
  tripletSpecificBaseQualityRatios[triplet,"max.alt"] = max.alt
  tripletSpecificBaseQualityRatios[triplet,"median.alt"] = median.alt
  tripletSpecificBaseQualityRatios[triplet,"n.alt"] = length(scores.alt)
  tripletSpecificBaseQualityRatios[triplet,"existingBaseScores.alt"] = paste0(names(table(scores.alt)), collapse = ",")
  
  tripletSpecificBaseQualityRatios[triplet,"ratio.mean"] = mean.alt / mean.ref
  tripletSpecificBaseQualityRatios[triplet,"ratio.median"] = median.alt / median.ref
  
  tripletSpecificBaseQualityRatios[triplet,"diff.mean"] = mean.alt - mean.ref
  tripletSpecificBaseQualityRatios[triplet,"diff.median"] = median.alt - median.ref  
  
  if(length(scores.alt)>1 & length(scores.ref)>1) {
    if (sum(abs(scores.alt - scores.ref)) == 0) {
      tripletSpecificBaseQualityRatios[triplet,"tTest_pval"] = 1
    } else {
      tripletSpecificBaseQualityRatios[triplet,"tTest_pval"] = t.test(scores.alt, scores.ref)$p.value
    }
  } else {
    if(length(scores.alt)>1 & length(scores.ref)==1) {
      tripletSpecificBaseQualityRatios[triplet,"tTest_pval"] = t.test(scores.alt, mu = scores.ref)$p.value
    } else if(length(scores.alt)==1 & length(scores.ref)>1) {
      tripletSpecificBaseQualityRatios[triplet,"tTest_pval"] = t.test(scores.ref, mu = scores.alt)$p.value
    } else if(length(scores.alt)==1 & length(scores.ref)==1) {
      tripletSpecificBaseQualityRatios[triplet,"tTest_pval"] = 1
    } else {
      tripletSpecificBaseQualityRatios[triplet,"tTest_pval"] = NA
      stop(paste0("something is strange here. We did not find a single value for REF and/or ALT for triplet",triplet,". Please check the base score files!"))
    }
  }
}

# adjust tTest pVals for multiple testing
tripletSpecificBaseQualityRatios$tTest_pval_adjusted = p.adjust(tripletSpecificBaseQualityRatios$tTest_pval, method = "fdr")


if (combineReverseComplement) {
  transitions = data.frame( c(rep("C",3),rep("T",3)), c("A","G","T","A","C","G"), stringsAsFactors = F)
} else {
  transitions = data.frame( c(rep("A",3),rep("C",3),rep("G",3),rep("T",3)), c("C","G","T","A","G","T","A","C","T","A","C","G"), stringsAsFactors = F)  
}
colnames(transitions) = c("FROM", "TO")



transition = function(value, start_point, end_point) {
  result = (value-start_point)/(end_point - start_point)
  result = max(min(1.0, result), 0)
  return(result)
}

# prepare color gradients which will be used in the plotting part downstream
col.ampel.red = colorRampPalette(c("white","white","red"), interpolate = "linear")(25+25)
col.ampel.green = colorRampPalette(c("green","white","white"), interpolate = "linear")(25+25)
col.ampel.greenWhiteRed = colorRampPalette(c("green","white","red"), interpolate = "linear")(40)


createAmpelColorFunction = function(col.ampel.type1, col.ampel.type2, intermediate.default, max.default) {
  # This function is a function generator. the result will be a function which can be used to get color gradients between 2 colors
  colorFunction = eval(parse(text = paste0("function(value, min=0, intermediate=",intermediate.default,", max=",max.default,") {
    color = ifelse(
      value<=intermediate,
      ",col.ampel.type1,"[round(transition(value, min, intermediate)*25+1,0)],
      ",col.ampel.type2,"[round((transition(value, intermediate, max)*25)+25,0)]
    )
    return (color)    
  }")))
  return(colorFunction)
}

getAmpelColor.red = createAmpelColorFunction(col.ampel.type1="col.ampel.red", col.ampel.type2="col.ampel.red", intermediate.default=100, max.default=1000)
getAmpelColor.green = createAmpelColorFunction(col.ampel.type1="col.ampel.green", col.ampel.type2="col.ampel.green", intermediate.default=100,max.default=1000)
getAmpelColor.greenWhiteRed = createAmpelColorFunction(col.ampel.type1="col.ampel.green", col.ampel.type2="col.ampel.red", intermediate.default=20,max.default=40)


## Add an alpha value to a colour
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) { rgb(x[1], x[2], x[3], alpha=alpha) }
  )  
}



######################################################
# PLOTTING PART STARTS HERE...

badBQ = 20 # values below which base quality will be called "poor quality" (red color)
intermediateBQ=32 # expected intermediate BQ (white color)
goodBQ=36 # values aboove which base quality will be called "good quality" (green color)



pdf(file = opt$outFile , width = 13, height = 5)
  # set opacity.pValThreshold to a value <0.99 in order to use opacity; 
  # e.g. fade out triplets changing insignificantly between ALT vs. REF base qualities  
  opacity.pValThreshold = 0.99 
  colorSaturation.lowerFoldChangeThreshold = 1.10 # start red/green coloring at this fold change
  colorSaturation.upperFoldChangeThreshold = 2.0  # full red/green coloring from this fold change on
  DiffThreshold = 2.00 # print difference values with abs value of at least colorSaturation.lowerDiffThreshold


  plot(c(0, 290), c(0, 40), type= "n", xlab = "", ylab = "", asp=1, xaxt='n', yaxt='n')
  text(130, 55,  opt$descriptionForMainTitle, cex = 1.1)
  text(-5, 5, "A", cex = 0.6);text(-5, 15, "C", cex = 0.6);text(-5, 25, "G", cex = 0.6);text(-5, 35, "T", cex = 0.6)
  
  
  xleft.baseAfter = c(0,10,20,30); names(xleft.baseAfter) = c("A","C","G","T")
  xright.baseAfter = xleft.baseAfter+10; names(xright.baseAfter) = c("A","C","G","T")
  ybottom.BaseBefore = c(0,10,20,30); names(ybottom.BaseBefore) = c("A","C","G","T")
  ytop.BaseBefore = ybottom.BaseBefore+10; names(ytop.BaseBefore) = c("A","C","G","T")
  transition.x_offset = c(0,50,100,150,200,250); names(transition.x_offset) = possible_mutations
  
  for (index.transition in seq(nrow(transitions))) {
    from = transitions[index.transition, "FROM"]
    to = transitions[index.transition, "TO"]
    mutation = paste0(from, to, collapse = "")
    text(transition.x_offset[mutation]+5, -5, "A", cex = 0.6);text(transition.x_offset[mutation]+15, -5, "C", cex = 0.6)
    text(transition.x_offset[mutation]+25, -5, "G", cex = 0.6);text(transition.x_offset[mutation]+35, -5, "T", cex = 0.6)
    for (baseBefore in c("A","C","G","T")) {
      for (baseAfter in c("A","C","G","T")) {
        triplet = paste0(baseBefore,from,to,baseAfter, collapse = "")
        indices = which(tripletSpecificBaseQualityRatios$triplet == triplet)
        xleft = transition.x_offset[mutation] + xleft.baseAfter[baseAfter]
        xright = transition.x_offset[mutation] + xright.baseAfter[baseAfter]
        ybottom = ybottom.BaseBefore[baseBefore]
        ytop = ytop.BaseBefore[baseBefore]      
        if (length(indices) > 0) {
  
          ratio.mean = tripletSpecificBaseQualityRatios[triplet,"ratio.mean"]
          diff.mean = tripletSpecificBaseQualityRatios[triplet,"diff.mean"]
          LFC.ratio = log2(ratio.mean)
          pVal = tripletSpecificBaseQualityRatios[triplet,"tTest_pval_adjusted"]
  
          if (plotType == Plottype.Ratio.UniColor) {
            # Color Ratios
            color = ifelse(LFC.ratio>=0,
                           getAmpelColor.red(value = LFC.ratio, min = 0.00, intermediate = log2(colorSaturation.lowerFoldChangeThreshold), max = log2(colorSaturation.upperFoldChangeThreshold)),
                           getAmpelColor.green(value = LFC.ratio, min = -log2(colorSaturation.upperFoldChangeThreshold), intermediate = -log2(colorSaturation.lowerFoldChangeThreshold), max = 0.00))

            alpha = transition(value = -log10(pVal), start_point = 0, end_point=-log10(opacity.pValThreshold))
            color = add.alpha(col = color, alpha = alpha)
            rect(xleft, ybottom, xright, ytop, col = color, border = "transparent")
            
            # write out ratio values
            if (LFC.ratio >= log2(colorSaturation.lowerFoldChangeThreshold) || LFC.ratio <= -log2(colorSaturation.lowerFoldChangeThreshold)) {
              text(xleft+(xright-xleft)/2, ybottom+(ytop-ybottom)/2, round(ratio.mean, 2), cex = 0.5)
            }
            
          } else if (plotType == Plottype.Differences.BiColor) {
            # Bi-Color by mean BQs
            meanBQ.ref = tripletSpecificBaseQualityRatios[triplet,"mean.ref"]
            meanBQ.alt = tripletSpecificBaseQualityRatios[triplet,"mean.alt"]
            alpha = transition(value = -log10(pVal), start_point = 0, end_point=-log10(opacity.pValThreshold))
            color.ref = getAmpelColor.greenWhiteRed(value = meanBQ.ref, min = badBQ, intermediate = intermediateBQ, max = goodBQ)
            color.ref = add.alpha(col = color.ref, alpha = alpha)
            color.alt = getAmpelColor.greenWhiteRed(value = meanBQ.alt, min = badBQ, intermediate = intermediateBQ, max = goodBQ)
            color.alt = add.alpha(col = color.alt, alpha = alpha)

            rect(xleft, ybottom, xleft+(xright-xleft)/2, ytop, col = color.ref, border = "transparent")
            rect(xleft+(xright-xleft)/2, ybottom, xright, ytop, col = color.alt, border = "transparent")
            rect(xleft, ybottom, xright, ytop, col = "transparent", border = "black")
            
            # write out difference values
            if (abs(diff.mean) >= DiffThreshold ) {
              text(xleft+(xright-xleft)/2, ybottom+(ytop-ybottom)/2, round(diff.mean, 2), cex = 0.5)
            }            
          }
          text(transition.x_offset[mutation] + 20, 45, paste0(from,"->",to), cex = 0.4)


        } else {
          rect(xleft, ybottom, xright, ytop, col = "grey", border = "transparent")
          rect(xleft, ybottom, xright, ytop, col = "transparent", border = "black")
        }
      }
    }
    rect(transition.x_offset[mutation]-0.3, -0.3, transition.x_offset[mutation]+40.3, 40.3, col = "transparent", border = "black")
  }

  
    
  if (plotType == Plottype.Differences.BiColor) {
    text(125, -14,  paste0("Values in cells show differences between ALT vs. REF base qualities"), cex = 0.6)
    text(125, -17,  paste0("pVal threshold: ",opacity.pValThreshold,
                           ";    green(BQ:",badBQ,")-white(BQ:",intermediateBQ,")-red(BQ:",goodBQ,");    BQ Difference threshold: +/-",DiffThreshold,";    BaseQuality threshold: ",opt$baseScoreThreshold), cex = 0.6)    
    
  } else if (plotType == Plottype.Ratio.UniColor) {
    text(125, -14,  paste0("Values in cells show ratios between ALT vs. REF base qualities"), cex = 0.6)
    text(125, -17,  paste0("pVal threshold: ",opacity.pValThreshold,
                           ";    starting to color green/red: ",round(1/colorSaturation.lowerFoldChangeThreshold,2),"/",colorSaturation.lowerFoldChangeThreshold," (LFC: -/+",round(log2(colorSaturation.lowerFoldChangeThreshold),2),")",
                           ";    fully green/red colored: ",round(1/colorSaturation.upperFoldChangeThreshold,2),"/",colorSaturation.upperFoldChangeThreshold," (LFC: -/+",round(log2(colorSaturation.upperFoldChangeThreshold),2),")",
                           ";    BaseQuality threshold: ",opt$baseScoreThreshold), cex = 0.6)    
  }


  
dev.off()  


write.table(tripletSpecificBaseQualityRatios, file=sub(".pdf", "_values.txt", opt$outFile), row.names = F, col.names = T, quote = F, sep = "\t")


# get set of all available base quality values
availabeBQs.ref = names(table(unlist(apply(tripletSpecificBaseQualityRatios, 1, function(line) {
  as.integer(unlist(strsplit(as.character(line["existingBaseScores.ref"]), ",")))
}))))
write.table(availabeBQs.ref, file=paste0(dirname(opt$outFile),"/availableBaseQualities.ref.txt"), quote = F, sep = "\t", row.names = F, col.names = F)

availabeBQs.alt = names(table(unlist(apply(tripletSpecificBaseQualityRatios, 1, function(line) {
  as.integer(unlist(strsplit(as.character(line["existingBaseScores.alt"]), ",")))
}))))
write.table(availabeBQs.alt, file=paste0(dirname(opt$outFile),"/availableBaseQualities.alt.txt"), quote = F, sep = "\t", row.names = F, col.names = F)

