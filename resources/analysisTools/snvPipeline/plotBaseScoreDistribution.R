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

numberOfCores = 5
library(doParallel)
registerDoParallel(numberOfCores)



## get parameters
opt = getopt(matrix(c(
  'vcfInputFile', 'v', 1, "character",
  'refScores', 'r', 1, "character",
  'altScores', 'a', 1, "character",
  'baseScoreThreshold', 't', 2, "integer",
  'descriptionForMainTitle', 'd', 2, "character",
  'outFile', 'o', 1, "character"
),ncol=4,byrow=TRUE));

 
checkForMissingParameter = function(parameter, errorText, exitCode=1) {
  if (is.null(opt[[parameter]])){
    cat(paste0(errorText,"\n")) 
    q(exitCode)      # quit, status unequal 0 means error
  }
}


checkForMissingParameter("vcfInputFile", "Please specify the file that contains the SNVs for which the base score distribution plot shall be created.", 1)
checkForMissingParameter("refScores", "Please specify the file containing reference allele base scores.", 1)
checkForMissingParameter("altScores", "Please specify the file containing alternative allele base scores.", 1)
checkForMissingParameter("outFile", "Please specify the output pdf file.", 1)


options(stringsAsFactors = FALSE)

locations = try( read.table(pipe(paste0("cut -f 1-2 ",opt$vcfInputFile)), comment.char = '', sep = "\t", header = T) )
  colnames(locations) = c("CHROM", "POS")
refScores = read.table(opt$refScores)
  colnames(refScores) = c("CHROM", "POS", "Base Score")
altScores = read.table(opt$altScores)
  colnames(altScores) = c("CHROM", "POS", "Base Score")


getWantedScores = function(scores, querySNVs, type) {
# return data frame with number of different BaseScore/VAF combinations for SNVs from scores which are also represented in querySNVs
  x=within(scores, PosString <- paste(CHROM,POS,sep='-'))[,"PosString"]
  y=within(querySNVs, PosString <- paste(CHROM,POS,sep='-'))[,"PosString"]
  indices = match(y,x)[!is.na(match(y,x))]
  wantedScores = scores[indices,]
  
  VAF = unlist(strsplit(wantedScores$`Base Score`, ';'))[c(F,T)]
  wantedScores$VAFclass = round(as.numeric(VAF) * 100/10,0)/10
  wantedScores$`Base Score` = unlist(strsplit(wantedScores$`Base Score`, ';'))[c(T,F)]
  wantedScores$`BaseScores` = strsplit(wantedScores$`Base Score`, ',')
  wantedScores = wantedScores[,-1:-3]

  indices = split(seq_len(nrow(wantedScores)), ceiling(seq_len(nrow(wantedScores))/ceiling(nrow(wantedScores)/numberOfCores/10)))

  data_frame_list_foreach_blocks =  foreach(i=1:length(indices)) %dopar% {
    data = wantedScores[indices[[i]],]
    tmp = as.data.frame(foreach(j=1:nrow(data), .combine=rbind) %do% {
      if (length(data[[j,"BaseScores"]])>1) {as.data.frame(sapply(data[j,], unlist))}
      else {as.data.frame(data[j,])}
    })
    tmp$BaseScores = as.integer(tmp$BaseScores)
    data.table(tmp)
  }

  resultTable = as.data.frame(table(rbindlist(data_frame_list_foreach_blocks)))
  resultTable$type = type

  return (resultTable)
}


getWantedScores_unaggregated = function(scores, querySNVs, type) {
# return data frame with BaseScor (BQ) and VAF for SNVs from scores which are also represented in querySNVs
  x=within(scores, PosString <- paste(CHROM,POS,sep='-'))[,"PosString"]
  y=within(querySNVs, PosString <- paste(CHROM,POS,sep='-'))[,"PosString"]
  indices = match(y,x)[!is.na(match(y,x))]
  wantedScores = scores[indices,]
  
  VAF = unlist(strsplit(wantedScores$`Base Score`, ';'))[c(F,T)]
  wantedScores$VAFclass = round(as.numeric(VAF) * 100/10,0)/10
  wantedScores$`Base Score` = unlist(strsplit(wantedScores$`Base Score`, ';'))[c(T,F)]
  wantedScores$`BaseScores` = strsplit(wantedScores$`Base Score`, ',')
  wantedScores = wantedScores[,-1:-3]

  indices = split(seq_len(nrow(wantedScores)), ceiling(seq_len(nrow(wantedScores))/ceiling(nrow(wantedScores)/numberOfCores/10)))
  
  data_frame_list_foreach_blocks =  foreach(i=1:length(indices)) %dopar% {
    data = wantedScores[indices[[i]],]
    tmp = as.data.frame(foreach(j=1:nrow(data), .combine=rbind) %do% {
      if (length(data[[j,"BaseScores"]])>1) {as.data.frame(sapply(data[j,], unlist))}
      else {as.data.frame(data[j,])}
    })
    tmp$BaseScores = as.integer(tmp$BaseScores)
    data.table(tmp)
  }

  resultTable <- do.call("rbind",data_frame_list_foreach_blocks)
  resultTable$type = type
  
  return (resultTable)
}


#Extract Legend 
g_legend<-function(a.gplot){ 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  return(legend)
} 

wantedRefScores = getWantedScores(scores = refScores, querySNVs = locations, type = "REF")
wantedAltScores = getWantedScores(altScores, locations, type = "ALT")

scores = rbind(wantedRefScores, wantedAltScores)
scores$type = factor(scores$type, levels = c("REF", "ALT"))
scores$BaseScores = as.integer(levels(scores$BaseScores))[scores$BaseScores]

threshold = opt$baseScoreThreshold

mainTitle = "Base Score distribution"
if (! is.null(opt$descriptionForMainTitle)) {
  mainTitle = paste0(mainTitle, "\n", opt$descriptionForMainTitle)
}

n=11
cbPalette = rainbow(n, s = 1, v = 1, start = 0, end = max(1, n - 1)/n, alpha = 1)

scores$VAFclassOrdered = factor(scores$VAFclass, levels = c(rev(levels(scores$VAFclass))))

pdf(file = paste0(opt$outFile), 8.27, 11.7) #paper = "a4"

    if (length(unique(scores$BaseScores)) < 10) {
      BaseScoreOffset = 0.4
      width = 0.75 
    } else {
      BaseScoreOffset = 0.2
      width = 0.37
    }
    scores$BaseScoresSplitted = ifelse(scores$type == "REF",scores$BaseScores-BaseScoreOffset,scores$BaseScores+BaseScoreOffset)
    
    density = ggplot(scores, aes(BaseScoresSplitted, Freq )) + geom_bar(aes(fill = VAFclassOrdered), width = width, position = "stack", stat="identity", alpha = 0.9) + theme(panel.margin = unit(0.125, "lines"))
    density = density + scale_fill_manual(values=cbPalette, guide = guide_legend(title = "VAFclass\n[REF/ALT]"), limits = seq(1, 0, -0.1))
    legend = g_legend(density)
    density = density + guides(fill=FALSE)
    if (length(unique(scores$BaseScores)) < 10) {
      density = density + scale_x_continuous(breaks=unique(scores$BaseScores))
    } else {
      density = density + scale_x_continuous(breaks = seq(min(scores$BaseScores), max(scores$BaseScores), 2))
    }
    density = density + ggtitle(paste0("Absolute ", mainTitle))
    if (! is.null(threshold)) {
      density = density + geom_vline(xintercept = threshold, colour="red")
    }
    density = density + theme(plot.title = element_text(hjust = 0.5))
    
    # calculate frequencies for relative representation;
    # relative means means that each Freq will be plotted relative to sum(Frequencies)[BaseScore].
    standardType = "REF"
    summedScores = dcast(scores, formula = BaseScores ~ type, value.var="Freq", fun=sum)
    scores = merge(scores, summedScores[,c("BaseScores",standardType)], by="BaseScores")
    scores$FreqNormalized = scores$Freq / scores$REF
    scores[scores$REF==0,"FreqNormalized"] = 0.0
    scores = scores[order(scores$BaseScores, scores$type, scores$VAFclassOrdered, decreasing = T),]

    relativePlot = ggplot(scores, aes(BaseScoresSplitted, Freq )) + geom_bar(aes(fill = VAFclassOrdered), width = width, position = "fill", stat="identity", alpha = 0.9) + theme(panel.margin = unit(0.125, "lines"))
    relativePlot = relativePlot + scale_fill_manual(values=cbPalette, guide = guide_legend(title = "VAFclass\n[REF/ALT]"), limits = seq(1, 0, -0.1)) + guides(fill=FALSE)
    if (length(unique(scores$BaseScores)) < 10) {
      relativePlot = relativePlot + scale_x_continuous(breaks=unique(scores$BaseScores))
    } else {
      relativePlot = relativePlot + scale_x_continuous(breaks = seq(min(scores$BaseScores), max(scores$BaseScores), 2))
    }
    relativePlot = relativePlot + ggtitle(paste0("Relative ", mainTitle))
    if (! is.null(threshold)) {
      relativePlot = relativePlot + geom_vline(xintercept = threshold, colour="red")
    }
    relativePlot = relativePlot + theme(plot.title = element_text(hjust = 0.5))

    relativePlot_facet = ggplot(scores, aes(BaseScores, Freq )) + geom_bar(aes(fill = VAFclassOrdered), width = 0.75, position = "fill", stat="identity", alpha = 0.9) + theme(panel.margin = unit(0.125, "lines")) +
      facet_wrap(~type, nrow=1)
    relativePlot_facet = relativePlot_facet + scale_fill_manual(values=cbPalette, guide = guide_legend(title = "VAFclass\n[REF/ALT]"), limits = seq(1, 0, -0.1)) + guides(fill=FALSE)
    if (length(unique(scores$BaseScores)) < 10) {
      relativePlot_facet = relativePlot_facet + scale_x_continuous(breaks=unique(scores$BaseScores))
    } else {
      relativePlot_facet = relativePlot_facet + scale_x_continuous(breaks = seq(min(scores$BaseScores), max(scores$BaseScores), 3))
    }
    relativePlot_facet = relativePlot_facet + ggtitle(paste0("Relative ", mainTitle))
    if (! is.null(threshold)) {
      relativePlot_facet = relativePlot_facet + geom_vline(xintercept = threshold, colour="red")
    }
    relativePlot_facet = relativePlot_facet + theme(plot.title = element_text(hjust = 0.5))
    
    relativePlot_RefNormalized = ggplot(scores, aes(BaseScoresSplitted, FreqNormalized )) + geom_bar(aes(fill = VAFclassOrdered), width = width, position = "stack", stat="identity", alpha = 0.9) + theme(panel.margin = unit(0.125, "lines"))    
    relativePlot_RefNormalized = relativePlot_RefNormalized + scale_fill_manual(values=cbPalette, guide = guide_legend(title = "VAFclass\n[REF/ALT]"), limits = seq(1, 0, -0.1)) + guides(fill=FALSE)
    if (length(unique(scores$BaseScores)) < 10) {
      relativePlot_RefNormalized = relativePlot_RefNormalized + scale_x_continuous(breaks=unique(scores$BaseScores))
    } else {
      relativePlot_RefNormalized = relativePlot_RefNormalized + scale_x_continuous(breaks = seq(min(scores$BaseScores), max(scores$BaseScores), 2))
    }    
    relativePlot_RefNormalized = relativePlot_RefNormalized + ggtitle(paste0("REF-relative ", mainTitle))
    if (! is.null(threshold)) {
      relativePlot_RefNormalized = relativePlot_RefNormalized + geom_vline(xintercept = threshold, colour="red")
    }
    relativePlot_RefNormalized = relativePlot_RefNormalized + theme(plot.title = element_text(hjust = 0.5))
    
    relativePlot_RefNormalized_facet = ggplot(scores, aes(BaseScores, FreqNormalized )) + geom_bar(aes(fill = VAFclassOrdered), width = 0.75, position = "stack", stat="identity", alpha = 0.9) + theme(panel.margin = unit(0.125, "lines")) +
      facet_wrap(~type, nrow=1) 
    relativePlot_RefNormalized_facet = relativePlot_RefNormalized_facet + scale_fill_manual(values=cbPalette, guide = guide_legend(title = "VAFclass\n[REF/ALT]"), limits = seq(1, 0, -0.1)) + guides(fill=FALSE)
    if (length(unique(scores$BaseScores)) < 10) {
      relativePlot_RefNormalized_facet = relativePlot_RefNormalized_facet + scale_x_continuous(breaks=unique(scores$BaseScores))
    } else {
      relativePlot_RefNormalized_facet = relativePlot_RefNormalized_facet + scale_x_continuous(breaks = seq(min(scores$BaseScores), max(scores$BaseScores), 3))
    }    
    relativePlot_RefNormalized_facet = relativePlot_RefNormalized_facet + ggtitle(paste0("REF-relative ", mainTitle))
    if (! is.null(threshold)) {
      relativePlot_RefNormalized_facet = relativePlot_RefNormalized_facet + geom_vline(xintercept = threshold, colour="red")
    }      
    relativePlot_RefNormalized_facet = relativePlot_RefNormalized_facet + theme(plot.title = element_text(hjust = 0.5))
    
    if (length(unique(scores$BaseScores)) < 10) {
      mytheme <- gridExtra::ttheme_default(
        core = list(fg_params=list(cex = 0.7)),
        colhead = list(fg_params=list(cex = 0.7)),
        rowhead = list(fg_params=list(cex = 0.7)))    
      
      numbers=(aggregate(Freq ~ BaseScores * type, scores, sum))
      numbers = t(dcast(numbers, formula = `BaseScores` ~ type, value.var="Freq"))
      colnames(numbers) = numbers["BaseScores",]
      numbers = numbers[-1,]
      
      numbersTable = qplot(1:10, 1:10, geom = "blank") + 
        theme_bw() +
        theme(line = element_blank(), text = element_blank()) +
        annotation_custom(grob = tableGrob(numbers, theme = mytheme), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)    
      
      
      # grid.newpage()
      gl <- grid.layout(nrow=7, ncol=2, widths = unit(c(3,1), "null"), heights = unit(c(3,3,3,3,3,1,1), "null"))
      pushViewport(viewport(layout = gl))
      
      vp.1 <- viewport(layout.pos.col=1, layout.pos.row=1) 
      vp.2 <- viewport(layout.pos.col=2, layout.pos.row=c(1,5)) 
      vp.3 <- viewport(layout.pos.col=1, layout.pos.row=2)
      vp.4 <- viewport(layout.pos.col=1, layout.pos.row=3)
      vp.5 <- viewport(layout.pos.col=1, layout.pos.row=4) 
      vp.6 <- viewport(layout.pos.col=1, layout.pos.row=5) 
      vp.7 <- viewport(layout.pos.col=c(1,2), layout.pos.row=c(6,7) )
      
      print(density, vp = vp.1)
      pushViewport(vp.2);grid.draw(legend);popViewport()
      print(relativePlot_RefNormalized, vp = vp.3)
      print(relativePlot_RefNormalized_facet, vp = vp.4)
      print(relativePlot, vp = vp.5)
      print(relativePlot_facet, vp = vp.6)
      
      print(numbersTable, vp = vp.7)
      
    } else {
      
      # grid.newpage()
      gl <- grid.layout(nrow=5, ncol=2, widths = unit(c(6,1), "null"), heights = unit(c(1,1,1,1,1), "null"))
      pushViewport(viewport(layout = gl))
      
      vp.1 <- viewport(layout.pos.col=1, layout.pos.row=1) 
      vp.2 <- viewport(layout.pos.col=2, layout.pos.row=c(1,5)) 
      vp.3 <- viewport(layout.pos.col=1, layout.pos.row=2)
      vp.4 <- viewport(layout.pos.col=1, layout.pos.row=3)
      vp.5 <- viewport(layout.pos.col=1, layout.pos.row=4) 
      vp.6 <- viewport(layout.pos.col=1, layout.pos.row=5) 
      
      print(density, vp = vp.1)
      pushViewport(vp.2);grid.draw(legend);popViewport()
      print(relativePlot_RefNormalized, vp = vp.3)
      print(relativePlot_RefNormalized_facet, vp = vp.4)
      print(relativePlot, vp = vp.5)
      print(relativePlot_facet, vp = vp.6)      
    }
    
dev.off()