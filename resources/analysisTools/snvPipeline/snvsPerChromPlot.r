#!/usr/bin/env Rscript
#
# Copyright (c) 2018 German Cancer Research Center (DKFZ).
#
# Distributed under the GPL-2 License (http://www.gnu.de/documents/gpl-2.0.en.html).
#

# to call this script: Rscript --vanilla /home/jaegern/ICGC_NGS/TOOLS/VARIANT_DETECTION/IntermutationDistance/snvsPerChromPlot.r -i intermutation_distance_snv_sorted.txt -s MB101 -o MB101_perChromSnvs


library(getopt)

plotColors = c("CA"="blue","CG"="black","CT"="magenta","TA"="purple","TC"="yellow","TG"="green")
baseChanges <- c("CA", "CG", "CT", "TA", "TC","TG")


# get parameters
opt = getopt(matrix(c(
  'mutationDistanceFile', 'i', 1, "character",
  'chromLengthFile', 'l', 1, "character",
  'sampleName', 's', 1, "character",  
  'outFilePrefix', 'o', 1, "character"
  ),ncol=4,byrow=TRUE));


if (is.null(opt$mutationDistanceFile)){      # no intermutation distance File specified
  cat("Please specify a tab separated intermutation distance file"); 
  q(status=1);      # quit, status unequal 0 means error
}
if (is.null(opt$chromLengthFile)){      # no chromLengthFile specified
  cat("Please specify the chromosome length file"); 
  q(status=1);      # quit, status unequal 0 means error
}
if (is.null(opt$sampleName)){
  cat("Please specify a sample name, which will be used as header for the plot");
  q(status=1);
}
if (is.null(opt$outFilePrefix)){
  cat("Please specify an outfile prefix");
  q(status=1);
}



dat = read.delim(opt$mutationDistanceFile, header = TRUE, stringsAsFactors = F)
# ensure standard chromosome order
dat$chromosome <- factor(dat$chromosome, levels = paste0("chr",c(seq(1,22),"X","Y")))


chromLength = read.table(file = opt$chromLengthFile, header = F)
rownames(chromLength) = paste0("chr",chromLength$V1)
chromLength = chromLength[,2, drop= F]
chromLengthMB = round(chromLength/1000000)


#chromMutRateNormalized = sapply(levels(dat$chromosome), function(x)(length(which(dat$chromosome == x))/ chromLengthMB[x]))
#barplot(chromMutRateNormalized, las=2, ylab="somatic snvs per MB", cex.lab=1.5)


baseChangeCount = sapply(levels(dat$chromosome), function(x) {
  sapply(baseChanges, function(y) {
    length(which(dat$mutationType[which(dat$chromosome == x)] == y))/chromLengthMB[x,]
  })
})


pdf(opt$outFilePrefix, width=9, height=6)
  # bar plot per chrom, normalized by chromosome length
  par(xpd=T, mar=par()$mar+c(0,0,0,3))    
  barplot(baseChangeCount, col = plotColors, main=paste(opt$sampleName), las=2, ylab="somatic snvs per MB", cex.lab=1.5)
  legend(29, (length(which(dat$chromosome == "chrX"))/ chromLengthMB["chrX",]), baseChanges[6:1], pch=15, col = plotColors[baseChanges[6:1]], cex=1, bty="n") 
dev.off()



