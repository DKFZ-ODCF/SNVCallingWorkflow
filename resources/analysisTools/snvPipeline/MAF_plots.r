
# R -f /home/jaegern/rWorkspace/MAF_plots.r --no-save --no-restore --args "sample_genomewide_high_conf_snvs_alleleFrq"

# Comparing Groups via Kernal Density 
# The sm.density.compare( ) function in the sm package allows to superimpose the kernal density plots of two or more groups.

library(getopt)
cmdArgs = commandArgs(TRUE)
print(cmdArgs)
if (length(cmdArgs) < 4) print(paste("Incorrect number of arguments (3 expected): ",length(cmdArgs)))
infile = cmdArgs[1]
snvnum = as.numeric(cmdArgs[2])
outfile = cmdArgs[3]
pid = cmdArgs[4]
nr_in_dbSNP = as.numeric(cmdArgs[5])
#library(Cairo)

maf=read.table(infile,sep="\t",as.is=TRUE)

pdf(outfile, width=7, height=7)

if(nr_in_dbSNP > 3){

	dens_in_dbSNP <- density(x = maf[maf[,1]==1,3])
	dens_in_all <- density(x = maf[,3])
	perc_in_dbSNP <- round((100*nr_in_dbSNP/snvnum),2)

	maxy <- max(dens_in_all$y)
	if(max(dens_in_all$y) < max(dens_in_dbSNP$y)){
		maxy <- max(dens_in_dbSNP$y)
	}

	plot(dens_in_dbSNP, ylim = c(0,maxy), col = "red", cex.lab=1.5, cex.axis=1.5, lwd=3, xlab=paste("Mutant allele frequency from ", snvnum, " SNVs and ", nr_in_dbSNP, " in dbSNP", sep=""), xlim=c(0, 1.0),main="")
	lines(dens_in_all,col = "blue", lwd = 3)
	abline(v=0.5,col="black",lty=2)
	title(main=paste(pid, perc_in_dbSNP, "% in dbSNP", sep=" "), cex.main = 1.5)

	legend("topright", c("ALL","IN dbSNP"), fill=c("blue","red"), cex = 0.7, lwd=1, bty="n")
}else{
#	dens_in_dbSNP <- density(x = maf[maf[,1]==1,3])
    dens_in_all <- density(x = maf[,3])
    perc_in_dbSNP <- round((100*nr_in_dbSNP/snvnum),2)

    maxy <- max(dens_in_all$y)
#    if(max(dens_in_all$y) < max(dens_in_dbSNP$y)){
#        maxy <- max(dens_in_dbSNP$y)
#    }

    plot(dens_in_all, ylim = c(0,maxy), col = "red", cex.lab=1.5, cex.axis=1.5, lwd=3, xlab=paste("Mutant allele frequency from ", snvnum, " SNVs and ", nr_in_dbSNP,     " in dbSNP", sep=""), xlim=c(0, 1.0),main="")
#    lines(dens_in_all,col = "blue", lwd = 3)
    abline(v=0.5,col="black",lty=2)
    title(main=paste(pid, perc_in_dbSNP, "% in dbSNP", sep=" "), cex.main = 1.5)

    legend("topright", c("ALL","IN dbSNP"), fill=c("blue","red"), cex = 0.7, lwd=1, bty="n")
}

dev.off()
