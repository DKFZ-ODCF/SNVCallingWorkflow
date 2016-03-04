#!/usr/bin/env Rscript

library(getopt)
library(ggplot2)
library(gridExtra) # for tableGrob
library(grid) # for gpar
library(reshape2) # for melt


## get parameters

opt = getopt(matrix(c(
  'vcfInputFile', 'v', 1, "character",
  'refScores', 'r', 1, "character",
  'altScores', 'a', 1, "character",
  'baseScoreThreshold', 't', 2, "integer",
  'descriptionForMainTitle', 'd', 2, "character",
  'outFile', 'o', 1, "character"
),ncol=4,byrow=TRUE));


if (is.null(opt$vcfInputFile)){
  cat("Please specify the file that contains the SNVs for which the base score distribution plot shall be created.\n"); 
  q(status=1);      # quit, status unequal 0 means error
}
if (is.null(opt$refScores)){
  cat("Please specify the file containing reference allele base scores.\n"); 
  q(status=1);      # quit, status unequal 0 means error
}
if (is.null(opt$altScores)){      # no vcf file specified
  cat("PPlease specify the file containing alternative allele base scores.\n"); 
  q(status=1);      # quit, status unequal 0 means error
}
if (is.null(opt$outFile)){      # no vcf file specified
  cat("Please specify the output pdf file.\n"); 
  q(status=2);      # quit, status unequal 0 means error
}


locations = try( read.table(pipe(paste0("cut -f 1-2 ",opt$vcfInputFile)), comment.char = '', sep = "\t", stringsAsFactors = F, header = T) )
colnames(locations) = c("CHROM", "POS")

refScores = read.table(opt$refScores)
  colnames(refScores) = c("CHROM", "POS", "Base Score")
  # refScores$type="REF"
altScores = read.table(opt$altScores)
  colnames(altScores) = c("CHROM", "POS", "Base Score")
  # altScores$type="ALT"

# filter base score entries for SNVs in input vcf file
wantedRefScores = refScores
wantedRefScores = merge(refScores, locations, by = c("CHROM","POS"))  
wantedRefScores = as.data.frame(unlist(strsplit(wantedRefScores[,3], ',')), stringsAsFactors = F)
colnames(wantedRefScores) = "Base Score"
wantedRefScores$`Base Score` = as.integer(wantedRefScores$`Base Score`)
wantedRefScores$type="REF"

wantedAltScores = altScores
wantedAltScores = merge(altScores, locations, by = c("CHROM","POS"))  
wantedAltScores = as.data.frame(unlist(strsplit(wantedAltScores[,3], ',')), stringsAsFactors = F)
colnames(wantedAltScores) = "Base Score"
wantedAltScores$`Base Score` = as.integer(wantedAltScores$`Base Score`)
wantedAltScores$type="ALT"
  
scores = rbind(wantedRefScores, wantedAltScores)
scores$type = factor(scores$type, levels = c("REF", "ALT"))


threshold = opt$baseScoreThreshold
mainTitle = "Base Score distribution"
if (! is.null(opt$descriptionForMainTitle)) {
  mainTitle = paste0(mainTitle, "\n", opt$descriptionForMainTitle)
}

numbers=table(scores)
cbPalette <- c("#56B4E9", "#CC79A7", "#009E73", "#D55E00", "#F0E442", "#0072B2", "#D55E00")
pdf(file = paste0(opt$outFile), paper = "a4")

  # binned base quality scores?
  # if (nrow(numbers) < 10) {
    numbers.m = melt(numbers, id.vars='Base Score')  

    density = ggplot(numbers.m, aes(`Base Score`, value )) + geom_bar(aes(fill = type), width = 0.65, position = "dodge", stat="identity", alpha = 0.9)
    density = density + scale_fill_manual(values=cbPalette)
    density = density + scale_x_discrete(breaks=rownames(numbers))
    density = density + ggtitle(mainTitle)
    if (! is.null(threshold)) {
      density = density + geom_vline(xintercept = threshold, colour="red")
    }
    
    if (nrow(numbers) < 10) {
      mytheme <- gridExtra::ttheme_default(
        core = list(fg_params=list(cex = 0.7)),
        colhead = list(fg_params=list(cex = 0.7)),
        rowhead = list(fg_params=list(cex = 0.7)))    
      
      numbersTable = qplot(1:10, 1:10, geom = "blank") + 
        theme_bw() +
        theme(line = element_blank(), text = element_blank()) +
        annotation_custom(grob = tableGrob(t(numbers), theme = mytheme), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)    
      
      grid.arrange(density, numbersTable, nrow = 2, heights = c(3,1))      
    } else {
      density
    }
      

#   } else {
#     density = ggplot(scores, aes(`Base Score`, fill = type)) + geom_density(alpha = 0.45)
#     density = density + scale_fill_manual(values=cbPalette)
#     density = density + ggtitle(mainTitle)
#     if (! is.null(threshold)) {
#       density = density + geom_vline(xintercept = threshold, colour="red")
#     }
#     
#     density
#   }
  
dev.off()

