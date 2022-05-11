#!/usr/bin/env python
#
# Copyright (c) 2018 German Cancer Research Center (DKFZ).
#
# Distributed under the MIT License (https://opensource.org/licenses/MIT).
#


# python /home/jaegern/pyWorkspace/NGS_Read_Processing/src/filter_PEoverlap.py --inf=snvs_108031.vcf --alignmentFile=/icgc/lsdf/mb/analysis/medullo/adultMB/results_per_pid/108031/alignment/tumor_108031_merged.bam.rmdup.bam --outf=snvs_108031_PEoverlapFiltered.vcf
# more snvs_108031.vcf | python /home/jaegern/pyWorkspace/NGS_Read_Processing/src/filter_PEoverlap.py --alignmentFile=/icgc/lsdf/mb/analysis/medullo/adultMB/results_per_pid/108031/alignment/tumor_108031_merged.bam.rmdup.bam --outf=snvs_108031_PEoverlapFiltered_nonALT_FINAL.vcf


import platform
if platform.python_implementation() == "PyPy":
    import copysam as pysam
else:  # "CPython"
    import pysam

import sys, os
from vcfparser import *
import re

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
    if remove_base.lower() == REF.lower():
        if remove_is_reverse:
            # check if none of the 4 DP4 values are < 0 now, which can happen due to BAQ values instead of original base qualities, which are not part of the BAM file
            if DP4rr > 0: DP4rr -= 1
        else:
            if DP4rf > 0: DP4rf -= 1
    elif remove_base.lower() == ALT.lower():
        if remove_is_reverse:
            if DP4ar > 0: DP4ar -= 1
        else:
            if DP4af > 0: DP4af -= 1   
    return(DP4rf, DP4rr, DP4af, DP4ar)

class BoolCounter:
    """ A class for counter objects """ 

    def __init__(self, ref_base, alt_base):
        self.ref_base = ref_base
        self.alt_base = alt_base
        self.ref_count = 0
        self.alt_count = 0

    def update(self, current_base):
        if self.ref_base == current_base:
            self.ref_count += 1
        elif self.alt_base == current_base:
            self.alt_count += 1

    def update_nonREFnonALT(self, count):
        self.alt_count =+ count

    def counted(self):
        return [self.ref_count, self.alt_count]

# MAIN ANALYSIS PROCEDURE
def performAnalysis(args):
    global qualScoreOffset

    # http://stackoverflow.com/questions/881696/unbuffered-stdout-in-python-as-in-python-u-from-within-the-program
    unbuffered = os.fdopen(sys.stdout.fileno(), 'w', 0)
    sys.stdout = unbuffered

    if args.qualityScore == 'illumina': qualScoreOffset = 64
    elif args.qualityScore == 'phred': qualScoreOffset = 33

    #vcfInFile = open(args.inf, "r")
    #outFile = open(args.outf, "w")

    # Reference file for BAQ_recalcuation and local realignment
    reference_file = pysam.Fastafile(args.refFileName)

    mode = "r"
    multiple_iterators = False
    # Setting pysam read mode based on the file extension
    if args.alignmentFile.split(".")[-1] == "bam":
        mode += "b"
    elif args.alignmentFile.split(".")[-1] == "cram":
        mode += "c"
    samfile = pysam.Samfile(args.alignmentFile, mode)  # This should work for BAM file only (with random access).

    if args.altPosF != '':
        ALT_basePositions_file = args.altPosF

    if args.refPosF != '':
        REF_basePositions_file = args.refPosF

    if args.altBQF != '':
        ALT_baseQualities_file = args.altBQF

    if args.refBQF != '':
        REF_baseQualities_file = args.refBQF


    for line in sys.stdin:  #   vcfInFile
        if line[0] == "#":
            headers = list(line[1:].rstrip().split('\t'))
            fixed_headers = ["^CHROM$", "^POS$", "^REF$", "^ALT$", "^INFO$", ]
            if args.no_control:
                fixed_headers += ["^CONFIDENCE$", "^RECLASSIFICATION$", ]
            else:
                fixed_headers += ["^ANNOTATION_control", ]
            header_indices = get_header_indices(headers, args.configfile, fixed_headers)
            sys.stdout.write(line)
            continue       # skip header from analysis

        entries = line.strip().split('\t')
        parsed_line = LineParser(entries, header_indices)

        nonREFnonALTfwd=0
        nonREFnonALTrev=0
        ALTcount=0

        ALT_basePositions=[]
        REF_basePositions=[]
        REF_baseQualities=[]
        ALT_baseQualities=[]

        # how to treat multiallelic SNVs? Skipped in this current version...                
        if ((args.no_control and int(parsed_line["CONFIDENCE"]) > 7 and "somatic" in parsed_line["RECLASSIFICATION"]) or (not args.no_control and "somatic" in parsed_line["ANNOTATION_control"])) and len(parsed_line["ALT"]) == 1:
            # DP=13;AF1=0.5;AC1=1;DP4=2,3,3,4;MQ=37;FQ=75;PV4=1,1,1,1
            info_values = parsed_line["INFO"].split(';')
            for info_idx, info_value in enumerate(info_values):
                if info_value[:4] == "DP4=":
                    DP4_idx = info_idx
                    DP4 = map(int, info_value[4:].split(','))
                    DP4rf, DP4rr, DP4af, DP4ar = DP4
                    DP4_original = re.sub('DP4', 'DP4original', info_value) # Keeping a backup of original DP4
                    DP4_original_alt = DP4af + DP4ar                    
                    break

            chrom=parsed_line["CHROM"]
            pos=int(parsed_line["POS"])
            REF=parsed_line["REF"]
            ALT=parsed_line["ALT"]
            
            readNameHash={}
            readMateHash={} # Hash to store read and mate starting positions for duplicate marking
            readMateHash_qnameLocation={} # Hash to store the location of gname in the above hash list


            ACGTNacgtn1 = [0]*10
            ACGTNacgtn2 = [0]*10
            count_PE = BoolCounter(REF, ALT) # Starting the counter for the forward and reverse reads removed due to PE overlap detection  
            count_supple = BoolCounter(REF, ALT) # "" for supplementary reads, since flag_filter is added, entire supplementary detection can be removed in future versions
            count_mismatch = BoolCounter(REF, ALT) # " for mismatch report 
            count_nonREFnonALT = BoolCounter(REF, ALT) # " to count the non-ref and non-alt base at POS

            # To match pysam and mpileup counts, a reference file is added. Given the reference file, Pysam by default computes BAQ (compute_baq).
            for pileupcolumn in samfile.pileup(chrom, (pos-1), pos, flag_filter=3844, redo_baq=True, ignore_overlaps=False, multiple_iterators=multiple_iterators):
                if pileupcolumn.pos == (pos-1):                	
                    #print 'coverage at base %s = %s' % (pileupcolumn.pos , pileupcolumn.nsegments)                    
                    for pileupread in pileupcolumn.pileups:                    	
                        if pileupread.is_del:                            
                            # 31 May 2016 JB: deletion at the pileup position                            
                            continue                        
                        baseScore = transformQualStr(pileupread.alignment.qual[pileupread.query_position])[0]
                        readpos = pileupread.query_position
                        if pileupread.alignment.seq[pileupread.query_position].lower()  == ALT.lower():                            
                            if args.altBQF != '':
                                ALT_baseQualities.append(baseScore)
                            if args.altPosF != '':
                                if pileupread.alignment.is_reverse:
                                    readlength = len(pileupread.alignment.seq)
                                    readpos = (readlength - readpos)
                                ALT_basePositions.append(readpos)
                        if pileupread.alignment.seq[pileupread.query_position].lower()  == REF.lower():
                            if args.refBQF != '':
                                REF_baseQualities.append(baseScore)
                            if args.refPosF != '':
                                if pileupread.alignment.is_reverse:
                                    readlength = len(pileupread.alignment.seq)
                                    readpos = (readlength - readpos)
                                REF_basePositions.append(readpos)

                        if pileupread.alignment.mapq >= args.mapq:
                            # http://wwwfgu.anat.ox.ac.uk/~andreas/documentation/samtools/api.html   USE qqual
                            try:
                                if transformQualStr(pileupread.alignment.qual[pileupread.query_position])[0] >= args.baseq:                                    
                                    # check if we consider this read as a proper read in terms of number of mismatches
                                    if args.allowedNumberOfMismatches > -1:
                                        numberOfMismatches = None
                                        for tag in pileupread.alignment.tags:
                                            if tag[0] == "NM":
                                                numberOfMismatches = tag[1]
                                                break 
                                            else:
                                                continue                                    
                                        
                                        if numberOfMismatches is not None:
                                            if numberOfMismatches > args.allowedNumberOfMismatches:
                                                remove_base = pileupread.alignment.seq[pileupread.query_position]
                                                remove_is_reverse = pileupread.alignment.is_reverse
                                                count_mismatch.update(remove_base)
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
                                        base = pileupread.alignment.seq[pileupread.query_position].lower()
                                        ACGTNacgtn_index = getIndexACGTNacgtn(is_reverse, is_read1, base)
                                        if(ACGTNacgtn_index[0] == "plus"):
                                            ACGTNacgtn1[ACGTNacgtn_index[1]] += 1
                                        else:
                                            ACGTNacgtn2[ACGTNacgtn_index[1]] += 1
 
                                        #if transformQualStr(pileupread.alignment.qual[pileupread.query_position])[0] >= args.baseq:        # DEBUG July 23 2012: BROAD BAM problem due to pileupread.alignment.qqual being shorter sometimes than pileupread.alignment.qual
                                        if(pileupread.alignment.query_name in readNameHash):
                                            #print pileupread.alignment.query_name
                                            old_qual = readNameHash[pileupread.alignment.query_name][0]
                                            old_base = readNameHash[pileupread.alignment.query_name][1]
                                            old_is_reverse = readNameHash[pileupread.alignment.query_name][2]
                                            old_read_mate_tuple = readNameHash[pileupread.alignment.query_name][3]
                                            current_qual = transformQualStr(pileupread.alignment.qual[pileupread.query_position])[0]
                                            current_base = pileupread.alignment.seq[pileupread.query_position]
                                            current_is_reverse = pileupread.alignment.is_reverse
                                            current_read_mate_tuple = (pileupread.alignment.reference_id, pileupread.alignment.reference_start, pileupread.alignment.reference_end, pileupread.alignment.next_reference_id, pileupread.alignment.next_reference_start)
                                            # if read name occurs twice for one variant, then due to overlapping PE reads, then subtract variant count from DP4 field
                                            # if old_base is not equal to new_base remove the one with the smaller base quality
                                            remove_base = None
                                            remove_is_reverse = None
                                            if(not(old_base == current_base)):
                                                if(old_qual <= current_qual):
                                                    remove_base = old_base
                                                    remove_is_reverse = old_is_reverse
                                                    remove_old = True
                                                else:
                                                    remove_base = current_base
                                                    remove_is_reverse = current_is_reverse
                                                    remove_old = False
                                            else:
                                                remove_base = current_base
                                                remove_is_reverse = current_is_reverse
                                                remove_old = False
                                            
                                            count_PE.update(remove_base)
                                            (DP4rf, DP4rr, DP4af, DP4ar) = decreaseDP4(remove_base, remove_is_reverse, REF, ALT, DP4rf, DP4rr, DP4af, DP4ar)                                                                                        
                                            # If current base is better, then removing the information about old mate
                                            # If current base is not good, then do nothing
                                            if remove_old:
                                            	old_location = readMateHash_qnameLocation[pileupread.alignment.query_name]
                                            	del readMateHash[old_read_mate_tuple][old_location]
                                            	read_mate_tuple_value = (current_qual, (current_base, current_is_reverse, pileupread.alignment.query_name, pileupread.alignment.is_supplementary))
                                            	if current_read_mate_tuple in readMateHash:
                                            		readMateHash[current_read_mate_tuple].append(read_mate_tuple_value)
                                            	else:
                                            		readMateHash[current_read_mate_tuple] = []
                                            		readMateHash[current_read_mate_tuple].append(read_mate_tuple_value)

                                        else:
                                            # Store base quality, base, and read direction in readNameHash
                                            base_qual_score=transformQualStr(pileupread.alignment.qual[pileupread.query_position])[0]
                                            read_mate_tuple = (pileupread.alignment.reference_id, pileupread.alignment.reference_start, pileupread.alignment.reference_end, pileupread.alignment.next_reference_id, pileupread.alignment.next_reference_start)
                                            read_mate_tuple_value = (base_qual_score, (pileupread.alignment.seq[pileupread.query_position], pileupread.alignment.is_reverse, pileupread.alignment.query_name, pileupread.alignment.is_supplementary))

                                            readNameHash[pileupread.alignment.query_name] = [base_qual_score, pileupread.alignment.seq[pileupread.query_position], pileupread.alignment.is_reverse, read_mate_tuple]

                                            if read_mate_tuple in readMateHash:
                                            	readMateHash[read_mate_tuple].append(read_mate_tuple_value)
                                            else:
                                            	readMateHash[read_mate_tuple] = []
                                            	readMateHash[read_mate_tuple].append(read_mate_tuple_value)

                                            readMateHash_qnameLocation[pileupread.alignment.query_name] = len(readMateHash[read_mate_tuple]) - 1 # Location of the last pushed element in the array                                            

                            except IndexError:
                                "soft-clipped or trimmed base, not part of the high-qual alignemnt anyways, skip"

                            if transformQualStr(pileupread.alignment.qual[pileupread.query_position])[0] >= args.baseq:
                            
                                if pileupread.alignment.seq[pileupread.query_position] == ALT:
                                    ALTcount += 1                                
                                # samtools mpileup sometimes counts bases as variants which are neither REF nor ALT
                                if (pileupread.alignment.seq[pileupread.query_position] != REF) and (pileupread.alignment.seq[pileupread.query_position] != ALT):                                    
                                    if pileupread.alignment.is_reverse:                                    	
                                        nonREFnonALTrev += 1
                                        #if DP4ar > 0: DP4ar -= 1
                                    else:
                                        nonREFnonALTfwd += 1
                                        #if DP4af > 0: DP4af -= 1

                    if (len(REF_baseQualities) > 0):
                        VAF =  1.0 * len(ALT_baseQualities) / (len(REF_baseQualities) + len(ALT_baseQualities))
                    elif (len(ALT_baseQualities) > 0):
                        VAF = 1.0
                    else:
                        VAF = 0.0

                    if args.altBQF != '':
                        scoreString = ",".join([str(score) for score in ALT_baseQualities])
                        if scoreString != '':
                            scoreString = ";".join([scoreString , str(VAF)])
                            ALT_baseQualities_file.write("%s\t%s\t%s\n" % (chrom, pos, scoreString))

                    if args.refBQF != '':
                        scoreString = ",".join([str(score) for score in REF_baseQualities])
                        if scoreString != '':
                            scoreString = ";".join([scoreString , str(VAF)])
                            REF_baseQualities_file.write("%s\t%s\t%s\n" % (chrom, pos, scoreString))

                    if args.altPosF != '':
                        positionsString = ",".join([str(readpos) for readpos in ALT_basePositions])
                        if positionsString != '':
                            ALT_basePositions_file.write("%s\t%s\t%s\n" % (chrom, pos, positionsString))

                    if args.refPosF != '':
                        positionsString = ",".join([str(readpos) for readpos in REF_basePositions])
                        if positionsString != '':
                            REF_basePositions_file.write("%s\t%s\t%s\n" % (chrom, pos, positionsString))

                    break # only one pileup for a position

            # Calculating duplicates based on read-mate pair's start positions (chr id and start location)
            count_duplicate = BoolCounter(REF, ALT)

            for key in readMateHash:
            	value_length = len(readMateHash[key])            	            	
            	if value_length > 0:            		
            	    sorted_values = sorted(readMateHash[key], key=lambda x: x[0]) # Sorted based on base quality
                    sorted_values = sorted_values[:-1] # removing the read with highest quality, so it will be retained for count
                    for value in sorted_values:        # Removing everthing else                            
                        qual_value, decreaseInfo = value
                        remove_base = decreaseInfo[0]
                        remove_is_reverse = decreaseInfo[1]

                        count_duplicate.update(remove_base)
                        (DP4rf, DP4rr, DP4af, DP4ar) = decreaseDP4(remove_base, remove_is_reverse, REF, ALT, DP4rf, DP4rr, DP4af, DP4ar)

            if (DP4[2] + DP4[3]) > ALTcount:    # that the ALTcount is larger  happens often due to BAQ during samtools mpileup which doesn't change the base qual in the BAM file, but decreases base qual during calling                
                if DP4af >= nonREFnonALTfwd: DP4af -= nonREFnonALTfwd

                if DP4ar >= nonREFnonALTrev: DP4ar -= nonREFnonALTrev

                nonREFnonALT = nonREFnonALTrev + nonREFnonALTfwd
                count_nonREFnonALT.update_nonREFnonALT(nonREFnonALT)

            #Format: DP2 -> "reference(forward + reverse), alt(forward + reverse)"
            supple_dup_str = 'DP2sup=' + ','.join(map(str, count_supple.counted()))
            supple_dup_str += ';DP2dup=' + ','.join(map(str, count_duplicate.counted()))
            supple_dup_str += ';DP2pairEnd=' + ','.join(map(str, count_PE.counted()))
            supple_dup_str += ';DP2mis=' + ','.join(map(str, count_mismatch.counted()))
            supple_dup_str += ';DP2nonREFnonALT=' + ','.join(map(str, count_nonREFnonALT.counted()))

            ACGTNacgtn1_string = "ACGTNacgtnPLUS="+",".join([str(i) for i in ACGTNacgtn1])
            ACGTNacgtn2_string = "ACGTNacgtnMINUS="+",".join([str(i) for i in ACGTNacgtn2])

            info_values[DP4_idx] = "DP4=" + str(DP4rf)+ "," + str(DP4rr)+ "," + str(DP4af)+ "," + str(DP4ar)            
            info_values.append(ACGTNacgtn1_string)
            info_values.append(ACGTNacgtn2_string)
            info_values.append(DP4_original)
            info_values.append(supple_dup_str)            

            entries[header_indices["INFO"]] = ';'.join(info_values)
            
            sys.stdout.write('\t'.join(entries) + '\n')
        else:
            sys.stdout.write(line)   # write germline and somatic-multiallelic SNVs as is

        del ALT_basePositions
        del REF_basePositions
        del REF_baseQualities
        del ALT_baseQualities
    samfile.close()

    if args.altBQF is not None:
        ALT_baseQualities_file.close()

    if args.refBQF is not None:
        REF_baseQualities_file.close()

    if args.altPosF is not None:
        ALT_basePositions_file.close()

    if args.refPosF is not None:
        REF_basePositions_file.close()

    #vcfInFile.close()
    #outFile.close()
    
if __name__ == '__main__':
    #print "Starting program...\n" 
    import argparse
    parser = argparse.ArgumentParser()
    #parser.add_option('--inf',action='store',type='string',dest='inf',help='Specify the name of the input vcf file containing all snvs (germline and somatic)',default='')
    parser.add_argument('--alignmentFile',dest='alignmentFile',help='Specify the name of the BAM file containing bwa alignments, has to be the BAM file that was used to call the variants in the input vcf file - REQUIRED', required=True)
    #parser.add_option('--outf',action='store',type='string',dest='outf',help='Specify the name of the output file, which will have same format as input vcf but with PE overlap filtered DP4 values if somatic and if snvs in PE overlap region',default='')
    parser.add_argument('--referenceFile', dest='refFileName', help='Specify the full path of reference genome file in fasta format, with index in the same directory')
    parser.add_argument('--mapq',type=int,dest='mapq',help='Specify the minimum mapping quality of bwa used for mpileup as parameter -q (default: 30 )',default=30)
    parser.add_argument('--baseq',type=int,dest='baseq',help='Specify the minimum base quality scores used for mpileup as parameter -Q (default: 13)',default=13)
    parser.add_argument('--qualityScore',dest='qualityScore',help='Specify whether the per base  quality score is given in phred or illumina format (default is Illumina score: ASCII offset of 64, while PHRED scores have an ASCII offset of 33)',default='phred')
    parser.add_argument('--maxNumberOfMismatchesInRead',type=int,dest='allowedNumberOfMismatches',help='Specify the number of mismatches that are allowed per read in order to consider this read. Value of -1 (default) turns this filter off.',default=-1)
    parser.add_argument('--altBaseQualFile',nargs="?",type=argparse.FileType('w'),dest='altBQF',help='Specify the name of the output file for alternative allele base qualities.',default=None)
    parser.add_argument('--refBaseQualFile',nargs="?",type=argparse.FileType('w'),dest='refBQF',help='Specify the name of the output file for reference allele base qualities.',default=None)
    parser.add_argument('--altBasePositionsFile',nargs="?",type=argparse.FileType('w'),dest='altPosF',help='Specify the name of the output file for position within the reads for alternative bases.',default=None)
    parser.add_argument('--refBasePositionsFile',nargs="?",type=argparse.FileType('w'),dest='refPosF',help='Specify the name of the output file for position within the reads for reference bases.',default=None)
    parser.add_argument('--nocontrol',action='store_true',dest='no_control',help='Specify the workflow is run without control.',default=False)
    parser.add_argument('--configFile',nargs="?",type=argparse.FileType('r'),dest='configfile',help='Specify the config file which contains header names.',default=None)


    args = parser.parse_args()
    performAnalysis(args)
    #print "\nProgram successfully terminating...."

