#!/usr/bin/env Rscript
#
# Copyright (c) 2018 German Cancer Research Center (DKFZ).
#
# Distributed under the GPL-2/GPL-3 License (http://www.gnu.de/documents/gpl-2.0.en.html, https://www.gnu.org/licenses/gpl-3.0.en.html).
#


cmdArgs = commandArgs(TRUE)
print(cmdArgs)
if (length(cmdArgs) < 3) print(paste("Incorrect number of arguments (3 expected): ",length(cmdArgs)))
infile = cmdArgs[1]
selection = cmdArgs[2]
outfile = cmdArgs[3]

library(methods)
library(ggplot2)
data <- read.table(infile, header=TRUE)
snvs <- data.frame(data)
sel <- which(snvs$type==selection)
somatic_snvs <- snvs[sel,]
all_somsnvs <- aggregate(data=somatic_snvs, value ~ REF + ALT + preceeding + following + type, sum)
all_somsnvs$REFALT <- paste(all_somsnvs$REF, all_somsnvs$ALT, sep="%->%")

#png(outfile, width=1200, height=300)
pdf(outfile, width=12, height=3)
p <- ggplot(all_somsnvs, aes(x=following,y=preceeding, fill=value))
p + geom_tile() + facet_grid( . ~ REFALT, labeller = label_parsed) + theme(panel.background = element_blank()) + scale_fill_gradient(low="yellow", high="red")
dev.off()
