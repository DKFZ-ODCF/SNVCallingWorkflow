#!/usr/bin/env python
#
# Copyright (c) 2018 German Cancer Research Center (DKFZ).
#
# Distributed under the MIT License (https://opensource.org/licenses/MIT).
#
#
# sort -k 1,1V -k 2,2n snvs_139625.annotated.somatic.vcf > snvs_139625_sort.vcf
# python /home/jaegern/pyWorkspace/NGS_Read_Processing/src/mutationDistance.py --inf=snvs_MB1.annotated.somatic.sort.vcf --outf=intermutation_distance_snv_TEST.txt


def listToTabsep(listItems, sep='\t'):
    return sep.join(listItems)

# strips off extra whitespace
def striplist(l):
    return([x.strip() for x in l])         


# MAIN ANALYSIS PROCEDURE
def performAnalysis(options):
    InFile = open(options.inf, "r")
    OutFile = open(options.outf, "w")
    OutFile.write("chromosome" + '\t' + "position" + '\t' + "intermutationDistance"+ '\t' + "mutationType" + '\t' + "genomicRegion" + '\t' + "MAF" +'\n')  
    excludedChromosomeList = options.excludedChromosomes.split(",")

    upstreamPos=0
    
    for line in InFile:
       	if line[0]=="#": continue 
        lineSplit=line.split('\t')
        lineSplitPlain=striplist(lineSplit)
        chrom=lineSplitPlain[0]
	if chrom.find("chr") == (-1): chrom="chr"+chrom
        if chrom in excludedChromosomeList: continue        
        DP4=map(int, (lineSplitPlain[7].split(';DP4=')[1].split(';')[0].split(',')) )
        percentVariantTumor=round(float((DP4[2]+DP4[3]))/sum(DP4),2)
        
        if percentVariantTumor < options.alleleFreq: continue  
        
        #if (lineSplitPlain[3]+lineSplitPlain[4] == "CT"):     
            
        pos=int(lineSplitPlain[1])
        snvDistance=pos - upstreamPos
        
        if snvDistance < 0:       # should only happen when chromosome changes
            snvDistance = pos
            
        if snvDistance == 0:      # can only happen when more than one sample as input, ie. the same snv in two samples
            snvDistance = 1


        baseChange=lineSplitPlain[3][0] + lineSplitPlain[4][0]
        if baseChange == "AC": baseChange = "TG"
        if baseChange == "AG": baseChange = "TC"        
        if baseChange == "AT": baseChange = "TA"
        if baseChange == "GT": baseChange = "CA"
        if baseChange == "GC": baseChange = "CG"        
        if baseChange == "GA": baseChange = "CT"

        genomicRegion = lineSplitPlain[15].split(';')[0]   # take annovar category as is for exonic, intergenic, intronic, all other annovar categories will be summarized under broader terms
        
        if genomicRegion[0:3] == "UTR": genomicRegion = "near_gene"   # if UTRs or down-/up-stream
        if genomicRegion.find("stream") > (-1): genomicRegion = "near_gene"
        if genomicRegion == "splicing": genomicRegion = "exonic"
        if genomicRegion[0:5] == "ncRNA": genomicRegion = "ncRNA" 
            
        
        OutFile.write(chrom + '\t' + str(pos) + '\t' + str(snvDistance) + '\t' + baseChange + '\t' + genomicRegion + '\t' + str(percentVariantTumor) + '\n')   
        upstreamPos=pos
    
    InFile.close()
    OutFile.close()
        

if __name__ == '__main__':
    print "Starting program...\n" 
    import optparse
    parser = optparse.OptionParser()
    parser.add_option('--inf',action='store',type='string',dest='inf',help='Specify the name of the input vcf file',default='')
    parser.add_option('--outf',action='store',type='string',dest='outf',help='Specify the name of the output file',default='')
    parser.add_option('--excludedChromosomes',action='store',type='string',dest='excludedChromosomes',help='Specify in a comma-separated string which chromosomes should be excluded', default='')
    parser.add_option('--alleleFreq',action='store',type='float',dest='alleleFreq',help='Specify the minimum mutant allele frequency for snvs to be used for analysis (default: 0 )',default=0)

    (options,args) = parser.parse_args()
    if len(options.inf) < 1 or options.outf == '':
        print "Mandatory parameters missing or wrong. Program will terminate now."
        print "\nYour parameter settings:"
        print options        
        raise SystemExit
    performAnalysis(options)     
    print "\nProgram successfully terminating...."  
