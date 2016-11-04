#!/usr/bin/env Rscript

# to call this script: Rscript --vanilla /home/jaegern/ICGC_NGS/TOOLS/VARIANT_DETECTION/IntermutationDistance/snvsPerChromPlot.r -i intermutation_distance_snv_sorted.txt -s MB101 -o MB101_perChromSnvs


library(getopt)

# get parameters

opt = getopt(matrix(c(
  'mutationDistanceFile', 'i', 1, "character",
  'sampleName', 's', 1, "character",  
  'outFilePrefix', 'o', 1, "character"
  ),ncol=4,byrow=TRUE));


if (is.null(opt$mutationDistanceFile)){      # no intermutation distance File specified
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


# main 
#png(opt$outFilePrefix, width=600, height=400)
pdf(opt$outFilePrefix, width=9, height=6)
# bar plot per chrom, normalized by chromosome length
par(xpd=T, mar=par()$mar+c(0,0,0,3))    

dat = read.delim(opt$mutationDistanceFile, header = TRUE)



plotColors = c("CA"="blue","CG"="black","CT"="magenta","TA"="purple","TC"="yellow","TG"="green")


# reorder to standrad chromosome order
dat$chromosome <- factor(dat$chromosome, levels = levels(dat$chromosome)[c(1,12,16,17,18,19,20,21,22,2,3,4,5,6,7,8,9,10,11,13,14,15,23)])

#plot(dat$chromosome, las=2)

#chromLength = read.table("/home/jaegern/rWorkspace/chromosomeLength.txt", header = TRUE)
chromLength = c("chr1"=249250621, "chr10"=135534747, "chr11"=135006516, "chr12"=133851895, "chr13"=115169878, "chr14"=107349540, "chr15"=102531392, "chr16"=90354753, "chr17"=81195210, "chr18"=78077248, "chr19"=59128983, "chr2"=243199373, "chr20"=63025520, "chr21"=48129895, "chr22"=51304566, "chr3"=198022430, "chr4"=191154276, "chr5"=180915260, "chr6"=171115067, "chr7"=159138663, "chr8"=146364022, "chr9"=141213431, "chrM"=16571, "chrX"=155270560, "chrY"=59373566) 
chromLengthMB = round(chromLength/1000000)


#chromMutRateNormalized = sapply(levels(dat$chromosome), function(x)(length(which(dat$chromosome == x))/ chromLengthMB[x]))
#barplot(chromMutRateNormalized, las=2, ylab="somatic snvs per MB", cex.lab=1.5)

baseChanges <- c("CA", "CG", "CT", "TA", "TC","TG")

#baseChangeCount = sapply(baseChanges, function(x)(length(which(dat$mutationType[chr1] == x))))

baseChangeCount = sapply(levels(dat$chromosome), function(x)(sapply(baseChanges, function(y)(length(which(dat$mutationType[which(dat$chromosome == x)] == y))/chromLengthMB[x]))))


barplot(baseChangeCount, col = plotColors, main=paste(opt$sampleName), las=2, ylab="somatic snvs per MB", cex.lab=1.5)

legend(29, (length(which(dat$chromosome == "chrX"))/ chromLengthMB["chrX"]), baseChanges[6:1], pch=15, col = plotColors[baseChanges[6:1]], cex=1, bty="n") 


dev.off()



