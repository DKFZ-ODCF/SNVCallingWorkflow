#!/usr/bin/env Rscript

library(getopt)
library(ggplot2)

## get parameters

opt = getopt(matrix(c(
  'refScores', 'r', 1, "character",
  'altScores', 'a', 1, "character",
  'baseScoreThreshold', 't', 2, "integer",
  'descriptionForMainTitle', 'd', 2, "character",
  'outFile', 'o', 1, "character"
),ncol=4,byrow=TRUE));


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


# refScores = as.data.frame(runif(1000, 50.0, 75))
# altScores = as.data.frame(runif(1000, 30.0, 60))
# threshold=42

refScores = read.table(opt$refScores)
  colnames(refScores) = "Base Score"
  refScores$type="REF"
altScores = read.table(opt$altScores)
  colnames(altScores) = "Base Score"
  altScores$type="ALT"

scores = rbind(refScores, altScores)
scores$type = factor(scores$type, levels = c("REF", "ALT"))

threshold = opt$baseScoreThreshold
mainTitle = "Base Score distribution"
if (! is.null(opt$descriptionForMainTitle)) {
  mainTitle = paste0(mainTitle, "\n", opt$descriptionForMainTitle)
}

pdf(file = paste0(opt$outFile), paper = "a4")
  p = ggplot(scores, aes(`Base Score`, fill = type)) + geom_density(alpha = 0.2)
  
  p = p + ggtitle(mainTitle)
  if (! is.null(threshold)) {
    p = p + geom_vline(xintercept = threshold, colour="red")
  }
  p
dev.off()
