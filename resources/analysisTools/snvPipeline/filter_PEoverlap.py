#! /usr/bin/python

# python /home/jaegern/pyWorkspace/NGS_Read_Processing/src/filter_PEoverlap.py --inf=snvs_108031.vcf --alignmentFile=/icgc/lsdf/mb/analysis/medullo/adultMB/results_per_pid/108031/alignment/tumor_108031_merged.bam.rmdup.bam --outf=snvs_108031_PEoverlapFiltered.vcf
# more snvs_108031.vcf | python /home/jaegern/pyWorkspace/NGS_Read_Processing/src/filter_PEoverlap.py --alignmentFile=/icgc/lsdf/mb/analysis/medullo/adultMB/results_per_pid/108031/alignment/tumor_108031_merged.bam.rmdup.bam --outf=snvs_108031_PEoverlapFiltered_nonALT_FINAL.vcf


import pysam
import sys


# strips off extra whitespace
def striplist(l):
    return([x.strip() for x in l])         


def listToTabsep(listItems, sep='\t'):
    return sep.join(listItems)


def qualFromASCII(ch):
    return(ord(ch) - qualScoreOffset)


def transformQualStr(s):
    return map(qualFromASCII,s)

# Routine for getting index of ACGTNacgtn lists
def getIndexACGTNacgtn(is_reverse, is_read1, base):
    if (is_reverse):
        if(is_read1):
            if(base == "a"):
                return ["minus", 5]
            elif(base == "c"):
                return ["minus", 6]
            elif(base == "g"):
                return ["minus", 7]
            elif(base == "t"):
                return ["minus", 8]
            elif(base == "n"):
                return ["minus", 9]
        else:
            if(base == "a"):
                return ["plus", 5]
            elif(base == "c"):
                return ["plus", 6]
            elif(base == "g"):
                return ["plus", 7]
            elif(base == "t"):
                return ["plus", 8]
            elif(base == "n"):
                return ["plus", 9]
    else:
        if(is_read1):
            if(base == "a"):
                return ["plus", 0]
            elif(base == "c"):
                return ["plus", 1]
            elif(base == "g"):
                return ["plus", 2]
            elif(base == "t"):
                return ["plus", 3]
            elif(base == "n"):
                return ["plus", 4]
        else:
            if(base == "a"):
                return ["minus", 0]
            elif(base == "c"):
                return ["minus", 1]
            elif(base == "g"):
                return ["minus", 2]
            elif(base == "t"):
                return ["minus", 3]
            elif(base == "n"):
                return ["minus", 4]

def decreaseDP4(remove_base, remove_is_reverse, REF, ALT, DP4rf, DP4rr, DP4af, DP4ar):
    if remove_base == REF:
        if remove_is_reverse:
            # check if none of the 4 DP4 values are < 0 now, which can happen due to BAQ values instead of original base qualities, which are not part of the BAM file
            if DP4rr > 0: DP4rr -= 1
        else:
            if DP4rf > 0: DP4rf -= 1
    elif remove_base == ALT:
        if remove_is_reverse:
            if DP4ar > 0: DP4ar -= 1
        else:
            if DP4af > 0: DP4af -= 1   
    return(DP4rf, DP4rr, DP4af, DP4ar)


# MAIN ANALYSIS PROCEDURE
def performAnalysis(options):
    global qualScoreOffset
    if options.qualityScore == 'illumina': qualScoreOffset = 64
    elif options.qualityScore == 'phred': qualScoreOffset = 33
    
    #vcfInFile = open(options.inf, "r")
    #outFile = open(options.outf, "w")


    samfile = pysam.Samfile(options.alignmentFile, "rb" ) # This should work for BAM file only (with random access).

    if options.altBQF != '':
        ALT_baseQualities_file = open(options.altBQF, 'w')

    if options.refBQF != '':
        REF_baseQualities_file = open(options.refBQF, 'w')


    for line in sys.stdin:           #   vcfInFile
        if line[0]=="#": 
            sys.stdout.write(line)
            continue       # skip header from analysis

        lineSplit=line.split('\t')
        lineSplitPlain=striplist(lineSplit)
        
        nonREFnonALTfwd=0
        nonREFnonALTrev=0
        ALTcount=0

        REF_baseQualities=[]
        ALT_baseQualities=[]
        
        if (lineSplitPlain[12].find("somatic") >= 0)  and (len(lineSplitPlain[4]) == 1):         # how to treat multiallelic SNVs? Skipped in this current version...
            # DP=13;AF1=0.5;AC1=1;DP4=2,3,3,4;MQ=37;FQ=75;PV4=1,1,1,1
            DP4start=lineSplitPlain[7].split('DP4=')[0]                  
            DP4=map(int,lineSplitPlain[7].split('DP4=')[1].split(';')[0].split(','))
            DP4rf=DP4[0]
            DP4rr=DP4[1]
            DP4af=DP4[2]
            DP4ar=DP4[3]
            
            sep=';'
            DP4end=sep.join(lineSplitPlain[7].split('DP4=')[1].split(';')[1:len(lineSplitPlain[7].split('DP4=')[1].split(';'))])
            
            chrom=lineSplitPlain[0]        
            pos=int(lineSplitPlain[1])
            REF=lineSplitPlain[3]
            ALT=lineSplitPlain[4]
            
            readNameHash={}

            ACGTNacgtn1 = [0]*10
            ACGTNacgtn2 = [0]*10

            for pileupcolumn in samfile.pileup(chrom, (pos-1), pos):
                if pileupcolumn.pos == (pos-1):
                    #print 'coverage at base %s = %s' % (pileupcolumn.pos , pileupcolumn.n)
                    for pileupread in pileupcolumn.pileups:
                        #print '\tbase in read %s = %s' % (pileupread.alignment.qname, pileupread.alignment.seq[pileupread.qpos])
                        baseScore = transformQualStr(pileupread.alignment.qual[pileupread.qpos])[0]
                        if pileupread.alignment.seq[pileupread.qpos].lower()  == ALT.lower():
                            ALT_baseQualities.append(baseScore)
                        if pileupread.alignment.seq[pileupread.qpos].lower()  == REF.lower():
                            REF_baseQualities.append(baseScore)
                        if pileupread.alignment.mapq >= options.mapq:

                            # http://wwwfgu.anat.ox.ac.uk/~andreas/documentation/samtools/api.html   USE qqual
                            try:
                                if transformQualStr(pileupread.alignment.qqual[pileupread.qpos])[0] >= options.baseq:
                                    # check if we consider this read as a proper read in terms of number of mismatches
                                    if options.allowedNumberOfMismatches > -1:
                                        numberOfMismatches = None
                                        for tag in pileupread.alignment.tags:
                                            if tag[0] == "NM":
                                                numberOfMismatches = tag[1]
                                                break 
                                            else:
                                                continue                                    
                                        
                                        if numberOfMismatches is not None:
                                            if numberOfMismatches > options.allowedNumberOfMismatches:
                                                remove_base = pileupread.alignment.seq[pileupread.qpos]
                                                remove_is_reverse = pileupread.alignment.is_reverse
                                                (DP4rf, DP4rr, DP4af, DP4ar) = decreaseDP4(remove_base, remove_is_reverse, REF, ALT, DP4rf, DP4rr, DP4af, DP4ar)                                         
                                                # after decreasing the respective DP4 value, go directly to the next read
                                                # without remembering the current read
                                                # This will lead to an unknown read name when the paired read occurs at the same
                                                # position. As we have already discarded the current high-mismatch read, we do not
                                                # have to decrease DP4 values again, when the read partner occurs at the same SNV.
                                                # We also do not increase ANCGTNacgtn for the discarded read.
                                                continue 
                                                                            
                                    # Check if pileupread.alignment is proper pair
                                    if(pileupread.alignment.is_proper_pair):
                                        # count to ACGTNacgtn list
                                        is_reverse = pileupread.alignment.is_reverse
                                        is_read1 = pileupread.alignment.is_read1
                                        base = pileupread.alignment.seq[pileupread.qpos].lower()
                                        ACGTNacgtn_index = getIndexACGTNacgtn(is_reverse, is_read1, base)
                                        if(ACGTNacgtn_index[0] == "plus"):
                                            ACGTNacgtn1[ACGTNacgtn_index[1]] += 1
                                        else:
                                            ACGTNacgtn2[ACGTNacgtn_index[1]] += 1

                                        #if transformQualStr(pileupread.alignment.qual[pileupread.qpos])[0] >= options.baseq:        # DEBUG July 23 2012: BROAD BAM problem due to pileupread.alignment.qqual being shorter sometimes than pileupread.alignment.qual
                                        if(readNameHash.has_key(pileupread.alignment.qname)):
                                            old_qual = readNameHash[pileupread.alignment.qname][0]
                                            old_base = readNameHash[pileupread.alignment.qname][1]
                                            old_is_reverse = readNameHash[pileupread.alignment.qname][2]
                                            current_qual = transformQualStr(pileupread.alignment.qqual[pileupread.qpos])[0]
                                            current_base = pileupread.alignment.seq[pileupread.qpos]
                                            current_is_reverse = pileupread.alignment.is_reverse
                                            # if read name occurs twice for one variant, then due to overlapping PE reads, then subtract variant count from DP4 field
                                            # if old_base is not equal to new_base remove the one with the smaller base quality
                                            remove_base = None
                                            remove_is_reverse = None
                                            if(not(old_base == current_base)):
                                                if(old_qual <= current_qual):
                                                    remove_base = old_base
                                                    remove_is_reverse = old_is_reverse
                                                else:
                                                    remove_base = current_base
                                                    remove_is_reverse = current_is_reverse
                                            else:
                                                remove_base = current_base
                                                remove_is_reverse = current_is_reverse
                                                
                                            (DP4rf, DP4rr, DP4af, DP4ar) = decreaseDP4(remove_base, remove_is_reverse, REF, ALT, DP4rf, DP4rr, DP4af, DP4ar)
                                        else:
                                            # Store base quality, base, and read direction in readNameHash
                                            readNameHash[pileupread.alignment.qname] = [transformQualStr(pileupread.alignment.qqual[pileupread.qpos])[0], pileupread.alignment.seq[pileupread.qpos], pileupread.alignment.is_reverse]
                            except IndexError:
                                "soft-clipped or trimmed base, not part of the high-qual alignemnt anyways, skip"
                                
                            
                            if transformQualStr(pileupread.alignment.qual[pileupread.qpos])[0] >= options.baseq:
                            
                                if pileupread.alignment.seq[pileupread.qpos] == ALT:
                                    ALTcount += 1
                                
                                # samtools mpileup sometimes counts bases as variants which are neither REF nor ALT
                                if (pileupread.alignment.seq[pileupread.qpos] != REF) and (pileupread.alignment.seq[pileupread.qpos] != ALT):
                                    if pileupread.alignment.is_reverse:
                                        nonREFnonALTrev += 1
                                        #if DP4ar > 0: DP4ar -= 1
                                    else:
                                        nonREFnonALTfwd += 1
                                        #if DP4af > 0: DP4af -= 1
                    break # only one pileup for a position

                    if options.altBQF != '':
                        scoreString = ",".join([str(score) for score in ALT_baseQualities])
                        if scoreString != '':
                            ALT_baseQualities_file.write("%s\t%s\t%s\n" % (chrom, pos, scoreString))

                    if options.refBQF != '':
                        scoreString = ",".join([str(score) for score in REF_baseQualities])
                        if scoreString != '':
                            REF_baseQualities_file.write("%s\t%s\t%s\n" % (chrom, pos, scoreString))


            if (DP4[2] + DP4[3]) > ALTcount:    # that the ALTcount is larger  happens often due to BAQ during samtools mpileup which doesn't change the base qual in the BAM file, but decreases base qual during calling
                #print line
                #print ALTcount
                #print (DP4[2] + DP4[3])
                if DP4af >= nonREFnonALTfwd: DP4af -= nonREFnonALTfwd
                if DP4ar >= nonREFnonALTrev: DP4ar -= nonREFnonALTrev
            
            ACGTNacgtn1_string = "ACGTNacgtnPLUS="+",".join([str(i) for i in ACGTNacgtn1])
            ACGTNacgtn2_string = "ACGTNacgtnMINUS="+",".join([str(i) for i in ACGTNacgtn2])
            DP4filtered = DP4start + "DP4=" + str(DP4rf)+ "," + str(DP4rr)+ "," + str(DP4af)+ "," + str(DP4ar)+ ";" + DP4end+";"+ACGTNacgtn1_string+";"+ACGTNacgtn2_string
            
            #if int(lineSplitPlain[28]) > 7:         #  filter only used in testing phase
            sys.stdout.write(listToTabsep(lineSplitPlain[0:7]) +'\t'+ DP4filtered +'\t'+ listToTabsep(lineSplitPlain[8:len(lineSplitPlain)]) +'\n')

        else:
            sys.stdout.write(line)   # write germline and somatic-multiallelic SNVs as is
    samfile.close()

    if options.altBQF != '':
        ALT_baseQualities_file.close()

    if options.refBQF != '':
        REF_baseQualities_file.close()
    #vcfInFile.close()
    #outFile.close()
    
    
if __name__ == '__main__':
    #print "Starting program...\n" 
    import optparse
    parser = optparse.OptionParser()
    #parser.add_option('--inf',action='store',type='string',dest='inf',help='Specify the name of the input vcf file containing all snvs (germline and somatic)',default='')
    parser.add_option('--alignmentFile',action='store',type='string',dest='alignmentFile',help='Specify the name of the BAM file containing bwa alignments, has to be the BAM file that was used to call the variants in the input vcf file - REQUIRED',default='')
    #parser.add_option('--outf',action='store',type='string',dest='outf',help='Specify the name of the output file, which will have same format as input vcf but with PE overlap filtered DP4 values if somatic and if snvs in PE overlap region',default='')
    parser.add_option('--mapq',action='store',type='int',dest='mapq',help='Specify the minimum mapping quality of bwa used for mpileup as parameter -q (default: 30 )',default=30)
    parser.add_option('--baseq',action='store',type='int',dest='baseq',help='Specify the minimum base quality scores used for mpileup as parameter -Q (default: 13)',default=13)
    parser.add_option('--qualityScore',action='store',type='string',dest='qualityScore',help='Specify whether the per base  quality score is given in phred or illumina format (default is Illumina score: ASCII offset of 64, while PHRED scores have an ASCII offset of 33)',default='phred')
    parser.add_option('--maxNumberOfMismatchesInRead',action='store',type='int',dest='allowedNumberOfMismatches',help='Specify the number of mismatches that are allowed per read in order to consider this read. Value of -1 (default) turns this filter off.',default=-1)
    parser.add_option('--altBaseQualFile',action='store',type='string',dest='altBQF',help='Specify the name of the output file for alternative allele base qualities.',default='')
    parser.add_option('--refBaseQualFile',action='store',type='string',dest='refBQF',help='Specify the name of the output file for reference allele base qualities.',default='')
    
    (options,args) = parser.parse_args()
    if len(options.alignmentFile) < 1:
        print "Mandatory parameters missing or wrong. Program will terminate now."
        print "\nYour parameter settings:"
        print options        
        raise SystemExit
    performAnalysis(options)     
    #print "\nProgram successfully terminating...."  

