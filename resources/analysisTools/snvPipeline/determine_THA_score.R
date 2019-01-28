#!/usr/bin/env Rscript
#
# Copyright (c) 2018 German Cancer Research Center (DKFZ).
#
# Distributed under the GPL-2 License (http://www.gnu.de/documents/gpl-2.0.en.html).
#
# This script determins the THA artifact score for the given snvs vcf file.
# The score is calculated as the fraction of SNVs with a baseQbias (PV4) p-value of at most 0.05 among all high quality somatic SNVs.
# See Phabricator Task 597 and related tasks.
# Author: G. Warsow

library(getopt)
options(stringsAsFactors = FALSE)

## get parameters
opt = getopt(matrix(c(
  'vcfInputFile', 'i', 1, "character"
),ncol=4,byrow=TRUE));

pval = 0.05

if (is.null(opt$vcfInputFile)){
  cat("Please specify the file that contains the somatic SNVs (conf8-10).\n"); 
  q(status=1);      # quit, status unequal 0 means error
}

PV4 = try( read.table(pipe(paste0("cat ",opt$vcfInputFile, " | cut -f 1,2,8 | perl -ne '$_ =~ /^(.*?)\t(.*?)\t.*?VDB=(.+?);.*?RPB=(.+?);.*?DP4=(.+?);.*?PV4=(.+?),(.+?),(.+?),(.+?);/; print \"$1\t$2\t$3\t$4\t$5\t$6\t$7\t$8\t$9\n\";' ")), comment.char = '', sep = "\t", header = T), silent = T )
if (class(PV4) == "try-error") {    
  stop("Could not read vcf file.")
  quit(1)
} 
colnames(PV4) = c("CHROM", "POS", "VDB", "RPB" , "DP4", "strandbias", "baseQbias", "mapQbias", "TDbias")
PV4[,paste0("baseQbias_below",pval)] = PV4$baseQbias < pval
count.BaseQ_below0.05 = sum(PV4[,paste0("baseQbias_below",pval)])
count.total = nrow(PV4)
# RELATIVE (not absolute) THA score:
THA_score = count.BaseQ_below0.05/count.total

cat(paste0(THA_score,"\n"))
