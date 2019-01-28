#!/usr/bin/env python
#
# Copyright (c) 2018 German Cancer Research Center (DKFZ).
#
# Distributed under the MIT License (https://opensource.org/licenses/MIT).
#

#
#  n.jaeger@dkfz.de
#  29 Nov 2011
#

import re
import sys

def listToTabsep(listItems, sep='\t'):
    return sep.join(listItems)


# MAIN ANALYSIS PROCEDURE
def performAnalysis(options):
    # bedInFile = open(options.inf, "r")
    bedOutFile = open(options.outf, "w")


    for line in sys.stdin:
        if line[0] == '#':
            bedOutFile.write(line)
            continue

        lineSplit=line.split('\t')
        try:
            DP4=lineSplit[7].split(';DP4=')[1].split(';')[0].split(',')
        except IndexError:
            print "Error while splitting: "
            print line
            sys.exit("Exiting...")

        ### General filters for SNVs and INDELS
        # DP4=104,47,0,1; don't keep snvs with just one read support
        if (DP4[2] == '0') and (DP4[3] == '1'): continue
        elif (DP4[2] == '1') and (DP4[3] == '0'): continue

        if DP4.count('0') >= 3:                    # skip snvs with only one entry in the DP4 field
            continue
        
        
        #DP4str=DP4
        DP4=map(int, DP4)
        
        if sum(DP4) < 4: continue

        leftToSnv,rightToSnv=lineSplit[-1].split(',')
        
        ### INDEL-specific filters
        if lineSplit[7][0:5] == 'INDEL':
            if float(lineSplit[5]) < 30.0:
                continue
            if float(lineSplit[5]) < 70.0:
                if ((DP4[2] == 0) or (DP4[3] == 0)) and ((DP4[2]+DP4[3]) <= 3):
                    continue
            leftToSnv += lineSplit[3]
            rightToSnv = lineSplit[3] + rightToSnv
        ### SNVs-specific filters
        else:
            if float(lineSplit[5]) < 3.0:
                percentVariantTumor=round(float((DP4[2]+DP4[3]))/sum(DP4),2)
                if percentVariantTumor < options.minVaf:  # skip if less than 5% of reads support variant 
                    continue
                if ((DP4[2]+DP4[3]) < options.minVac):
                    continue
            
            if float(lineSplit[5]) < 20.0:
                if ((DP4[2] == 0) or (DP4[3] == 0)) and ((DP4[2]+DP4[3]) <= options.minVacPS):
                    continue

        ### Strand-bias filters for Indels and SNVs (for Indels, the REF sequence is added to the sequence context to test; maybe this is too conservative for deletions...)
        if (DP4[3] < 2):       # forward strand bias, i.e. less than two reads on reverse strand for variant
            pattern = "GG(.){2}GG" 
            if (re.search(pattern, leftToSnv) or (leftToSnv.find('GGC') > -1)): # filter GGNNGG motif and GGC base triplet from NAR publication "Sequence-specific error profile of Illumina sequencers"
                if (DP4[3] != 1) and (DP4[2] != 2):
                    continue
                
        if (DP4[2] < 2):       # reverse strand bias, , i.e. less than two reads on forward strand for variant
            pattern = "CC(.){2}CC"
            if (re.search(pattern, rightToSnv) or (rightToSnv.find('GCC') > -1)): 
                if (DP4[2] != 1) and (DP4[3] != 2):
                    continue


        bedOutFile.write(line)   


    # bedInFile.close()
    bedOutFile.close()
        

if __name__ == '__main__':
    print "Starting program...\n" 
    import optparse
    parser = optparse.OptionParser()
    # parser.add_option('--inf',action='store',type='string',dest='inf',help='Specify the name of the input bed file',default='')
    parser.add_option('--outf',action='store',type='string',dest='outf',help='Specify the name of the output bed file',default='')
    parser.add_option('--minVaf',action='store',type='float',dest='minVaf',help='Specify the minimum variant allele frequency for variants with a bcf score below 3.0 (default: 0.05 )',default=0.05)
    parser.add_option('--minVac',action='store',type='int',dest='minVac',help='Specify the minimum count of variant reads for variants with a bcf score below 3.0 (default: 5 )',default=5)
    parser.add_option('--minVacPS',action='store',type='int',dest='minVacPS',help='Specify the minimum count of variant reads for one strand if the other strand has 0 reads for variants with a bcf score below 20.0 (default: 3 )',default=3)
    (options,args) = parser.parse_args()
    if options.outf == '':
        print "Mandatory parameters missing or wrong. Program will terminate now."
        print "\nYour parameter settings:"
        print options        
        raise SystemExit
    performAnalysis(options) 
    print "\nProgram successfully terminating...."
