#!/usr/bin/env Rscript
#
# Copyright (c) 2018 German Cancer Research Center (DKFZ).
#
# Distributed under the GPL-2/GPL-3 License (http://www.gnu.de/documents/gpl-2.0.en.html, https://www.gnu.org/licenses/gpl-3.0.en.html).
#

### to call this script: Rscript --vanilla ~/org/Bioinfo/rainfallPlots/intermutationDistance_Coord_color.r -i snvFile.vcf -s MB1000 -o MB1000_intermutationDistance_snvs
### pay attention with the chrLengthFile; there is a default which is either set to hg19 or hs37d5 

library(getopt)

## get parameters

opt = getopt(matrix(c(
  'mutationFile', 'i', 1, "character",
  'sampleName', 's', 1, "character",  
  'outFilePrefix', 'o', 1, "character",
  'chrLengthFile', 'l', 1, "character",
  'chrArray', 'a', 1, "character",
  'chrPrefix', 'p', 2, "character",
  'chrSuffix', 'u', 2, "character"

  ),ncol=4,byrow=TRUE));

if (is.null(opt$mutationFile)){      # no intermutation distance File specified
  cat("Please specify a tab separated intermutation distance file"); 
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

if (is.null(opt$chrLengthFile)){
  opt$chrLengthFile = '/ibios/co02/reference/hg19/hg19_chrTotalLength.tsv'
}

c.array <- unlist(strsplit(opt$chrArray, ',', fixed = TRUE))

data <- read.table(opt$chrLengthFile, header=FALSE)
chrLength <- data.frame(data)
rownames(chrLength) <- chrLength$V1
xtotal <- sum(chrLength[c.array,2]/10)
xoffsets <- c(0)

headset <- read.delim(pipe(paste0("grep -v '^##' ",opt$mutationFile)), header = TRUE)
headset.classes = sapply(headset, class)
headset.classes[names(headset.classes) %in% c("REF", "ALT")] = "character"

dat <- read.delim(pipe(paste0("grep -v '^##' ",opt$mutationFile)), header = TRUE, colClasses = headset.classes)
dat <- dat[which(dat$REF != 'N'),]
if (nrow(dat) > 0) {
  dat$diff <- c(-1,diff(dat$POS))
}

complement <- c('A' = 'T', 'C' = 'G','G' = 'C', 'T' = 'A', 'N' = 'N')

sel <- which(dat$REF == 'C' | dat$REF == 'T')
if (length(sel) > 0) {
  dat[sel,'change'] <- paste0(dat[sel, 'REF'],dat[sel, 'ALT'])
}

sel <- which(dat$REF == 'A' | dat$REF == 'G')
if (length(sel) > 0) {
  dat[sel, 'change'] <- paste0(complement[dat[sel, 'REF']],complement[dat[sel, 'ALT']])
}


##chrX <- which(dat$chromosome == "X")
##chrXmean <- round(mean(dat$intermutationDistance[chrX])/1000)  # in kb

##chrAllbutX <- which(dat$chromosome != "X")
##chrAllmean <- round(mean(dat$intermutationDistance[chrAllbutX])/1000) # in kb

### main 
#png(paste(opt$outFilePrefix,".png",sep = ""), width=1000, height=400)
pdf(opt$outFilePrefix, width=20, height=8)

plotColors = c("CA"="blue","CG"="black","CT"="red","TA"="purple","TC"="orange","TG"="green")
baseChanges <- c("CA", "CG", "CT", "TA", "TC","TG")

c.name <- paste0(opt$chrPrefix, c.array[1], opt$chrSuffix)  

sel <- which(dat$diff > 0)
if(length(sel)){
	maxY <- max(log10(dat$diff[sel]))

	sel <- which(dat$X.CHROM == c.name & dat$diff > 0)

	plot(0, type='n', main=opt$sampleName, xlab="",ylab="intermutation distance (bp)",cex.lab=1.5, xlim=c(0.02,1.035), ylim=c(-0.05,maxY), axes=FALSE)

	#xoffset <- chrLength[chrLength$V1==c.array[1],2]/10
	xoffset <- 0


	# chromosomes <- c(2:22, 'X')

	for (chr in c.array) {
	  c.name <- paste0(opt$chrPrefix, chr, opt$chrSuffix)  

	  sel <- which(dat$X.CHROM == c.name & dat$diff > 0)
	  points((dat$POS[sel]/10+xoffset)/xtotal, log10(dat$diff[sel]), col = plotColors[dat$change[sel]], pch=16,cex=0.8)
	  
	  xoffset <- xoffset + chrLength[chrLength$V1==chr,2] / 10
	  xoffsets <- append(xoffsets, xoffset)
	  
	}


	abline(v=(xoffsets/xtotal), col="#000000")

	axis(2,at=c(0,2,4,6),format(c(1,100,10000,1000000), scientific=FALSE), cex.axis = 1.5)  
	box()
	sel <- which(dat$diff > 0)
	abline(h=(seq(0, max(log10(dat$diff[sel])), 1)), col="lightgray", lty="dotted", lwd=2)


	textXpos <- chrLength[c.array,2] / 20 + xoffsets[1:length(xoffsets) -1 ]
	text(textXpos/xtotal, 0, labels = paste("",c.array,sep=''), pos = 1, cex=0.8)

	legend("right", baseChanges, cex=1.3, fill=plotColors[baseChanges], bty="n")
}else{
	plot(0, type='n', main=opt$sampleName, xlab="",ylab="intermutation distance (bp)",cex.lab=1.5, xlim=c(0.02,1.035), ylim=c(-0.05,10000), axes=FALSE)
	text(0.5, 5000, labels = "There is no chromosome with more than one somatic SNV in the data", pos = 1, cex=2)
}
dev.off()
