#!/usr/bin/env python
#
# Copyright (c) 2018 German Cancer Research Center (DKFZ).
#
# Distributed under the MIT License (https://opensource.org/licenses/MIT).
#
#MS, 2015-01
#Adapted from MB,11.2014 and OE,11.2014
import sys
import optparse
import pysam
import numpy as np
import math
# import matplotlib
# matplotlib.use('Agg')
# import matplotlib.pyplot as plt
# import matplotlib.gridspec as gridspec
import re
from math import exp, expm1
from Bio.Seq import Seq
from scipy.stats import binom


def readSeqBiasMatrix(FN):
    File = open(FN, "r")
    error_matrix = {}
    mutation_count_matrix = {}
    for line in File:
        if(line[0] == "#"):
            mutation = line[1:].rstrip().split("\t")
            error_matrix[mutation[0]] = {}
            mutation_count_matrix[mutation[0]] = {}
        elif(line[0] == "\t"):
            header = line[1:].rstrip().split("\t")
        else:
            split_line = line.rstrip().split("\t")
            error_matrix[mutation[0]][split_line[0]]={}
            mutation_count_matrix[mutation[0]][split_line[0]]={}
            j=0    #index of mutation_after header
            for i in range(1,len(split_line)):
                match = re.search(r'(\d+)\/(\d+);(\d+)',split_line[i])
                error_matrix[mutation[0]][split_line[0]][header[j]] = [match.group(1),match.group(2)]
                mutation_count_matrix[mutation[0]][split_line[0]][header[j]] = match.group(3)
                j=j+1
                if(j == 4):
                    j=0
    return error_matrix, mutation_count_matrix



def calculateBiasMatrix(error_matrix, mutation_count_matrix,numReads,numMuts,biasPValThreshold,biasRatioThreshold):

    possible_mutations = ["CA", "CG", "CT", "TA", "TC", "TG"]

    bias_matrix = {}
    for mutation in possible_mutations:
        bias_matrix[mutation] = {}
        possible_bases = ["A", "C", "G", "T"]
        for base_before in possible_bases:
            bias_matrix[mutation][base_before] = {}
            for base_after in possible_bases:
                bias_matrix[mutation][base_before][base_after] = 0
                nReads=int(error_matrix[mutation][base_before][base_after][0])+int(error_matrix[mutation][base_before][base_after][1])
                if(nReads > 0):
                    sign=int(math.copysign(1,int(error_matrix[mutation][base_before][base_after][0])-int(error_matrix[mutation][base_before][base_after][1])))
                    fracPlus=float(error_matrix[mutation][base_before][base_after][0])/nReads
                    fracMinus=float(error_matrix[mutation][base_before][base_after][1])/nReads
                    minorReadCount = int(min(int(error_matrix[mutation][base_before][base_after][0]),int(error_matrix[mutation][base_before][base_after][1])))
                    p = binom.cdf(minorReadCount, nReads, 0.5)
                    nMuts=mutation_count_matrix[mutation][base_before][base_after]
#                    print str(mutation) + " in " + str(base_before) + str(mutation[0]) + str(base_after) + ": " + str(p) + " (" + str(minorReadCount) + ";" + str(nReads) + ")"
                    if(p < biasPValThreshold and (fracPlus > biasRatioMinimum or fracMinus > biasRatioMinimum)):
                        if(nReads > numReads and nMuts > numMuts):
                            if(fracPlus > biasRatioThreshold or fracMinus > biasRatioThreshold):
                                sign *= 2
                            bias_matrix[mutation][base_before][base_after] = sign

    return bias_matrix

def printErrorMatrix(error_matrix, output_filename):
    possible_mutations = ["CA", "CG", "CT", "TA", "TC", "TG"]
    bases = ["A", "C", "G", "T"]
    error_file=open(output_filename, "w")
    possible_mutations=sorted(error_matrix.keys())
    for mutation in possible_mutations:
        error_file.write("#"+mutation+"\n")
        error_file.write("\t".join([""]+bases)+"\n")
        for base in bases:
            error_list = [str(error_matrix[mutation][base][after]) for after in bases]
            error_file.write("\t".join([base]+error_list)+"\n")
    error_file.close()





def complement(base):
    if(base == "A"):
        return "T"
    elif(base == "C"):
        return "G"
    elif(base == "G"):
        return "C"
    elif(base == "T"):
        return "A"
    elif(base == "a"):
        return "t"
    elif(base == "c"):
        return "g"
    elif(base == "g"):
        return "c"
    elif(base == "t"):
        return "a"
    elif(base == "n"):
        return "n"
    elif(base == "N"):
        return "N"
    else:
        return base

def markSeqBiasInVcfFile(vcfFilename, referenceFilename, bias_matrixSeq, bias_matrixSeqing, vcfFlaggedFilename, maxNumOppositeReadsSequenceWeakBias, maxNumOppositeReadsSequencingWeakBias, maxNumOppositeReadsSequenceStrongBias, maxNumOppositeReadsSequencingStrongBias, ratioVcfSequence, ratioVcfSequencing):
    vcfFile = open(vcfFilename, "r")
    vcfFlagged=open(vcfFlaggedFilename, "w")
    reference = None
    if(not(referenceFilename == "NA")):
        reference = pysam.Fastafile(referenceFilename)

    possible_mutations = ["CA", "CG", "CT", "TA", "TC", "TG"]

    header = []
    for line in vcfFile:
        if line[0:2] == "##":
            vcfFlagged.write(line)
            continue
        elif line[0] == "#":
            header = line[1:].rstrip().split("\t")
            vcfFlagged.write(line.rstrip()+"\t"+"seqBiasPresent"+"\t"+"seqingBiasPresent"+"\n")
        #elif line.rstrip().split("\t")[header.index("ANNOTATION_control")] != "somatic":
        elif not "ACGTNacgtnPLUS" in line or not "ACGTNacgtnMINUS" in line:
            vcfFlagged.write(line.rstrip()+"\t"+"."+"\t"+"."+"\n")
            continue
        else:
            seqBiasPresent=0
            seqingBiasPresent=0
            possible_bases = ["A", "C", "G", "T", "N", "a", "c", "g", "t", "n"]
            posible_proper_bases = ["A", "C", "G", "T", "a", "c", "g", "t"]
            split_line = line.rstrip().split("\t")

            chrom = split_line[header.index("CHROM")]
            pos = int(split_line[header.index("POS")])
            context = ""
            if(referenceFilename == "NA"):
                context = split_line[header.index("SEQUENCE_CONTEXT")].split(",")[0][-1]+split_line[header.index("REF")]+split_line[header.index("SEQUENCE_CONTEXT")].split(",")[1][0]
            else:
                context = reference.fetch(chrom, pos-2, pos+1)

            current_mutation = split_line[header.index("REF")]+split_line[header.index("ALT")]
            base_before = context[0].upper()
            base_after = context[2].upper()

            # Do not test variants with non-ACGTN flanking bases; just write them out
            if base_before not in posible_proper_bases or base_after not in posible_proper_bases:
                vcfFlagged.write(line.rstrip()+"\t"+"."+"\t"+"."+"\n")
                continue

            info_list = [i.split("=") for i in split_line[header.index("INFO")].split(";")]

            # Get strand specific counts
            ACGTNacgtnPLUS = []
            ACGTNacgtnMINUS = []

            for element in info_list:
                if(element[0] == "ACGTNacgtnPLUS"):
                    ACGTNacgtnPLUS = [int(i) for i in element[1].split(",")]
                elif(element[0] == "ACGTNacgtnMINUS"):
                    ACGTNacgtnMINUS = [int(i) for i in element[1].split(",")]

            # Count number of alt bases
            read1_nr = ACGTNacgtnPLUS[possible_bases.index(current_mutation[1])]
            read1_r = ACGTNacgtnMINUS[possible_bases.index(current_mutation[1].lower())]
            read2_nr = ACGTNacgtnMINUS[possible_bases.index(current_mutation[1])]
            read2_r = ACGTNacgtnPLUS[possible_bases.index(current_mutation[1].lower())]



            PCR_plus = 0
            PCR_minus = 0

            SEQ_plus = 0
            SEQ_minus = 0

            # reverse=0

            try:
                mutation_index = possible_mutations.index(current_mutation)
                #REF base is C or T
                PCR_plus = read1_nr + read2_r
                PCR_minus = read2_nr + read1_r

                SEQ_plus = read1_nr + read2_nr
                SEQ_minus = read1_r + read2_r

                # reverse=0

            except ValueError:
                # REF base is A or G -> reverse complement
                current_mutation = complement(current_mutation[0])+complement(current_mutation[1])
                base_before_reverse_complement = complement(base_after)
                base_after_reverse_complement = complement(base_before)

                base_before = base_before_reverse_complement
                base_after = base_after_reverse_complement

                # TODO: is this correct or do we need to rev comp? Changed it!
                PCR_minus  = read1_nr + read2_r
                PCR_plus = read2_nr + read1_r

                SEQ_minus = read1_nr + read2_nr
                SEQ_plus = read1_r + read2_r

                # reverse=1


            #these totals are equal
            PCRTot=float(PCR_plus)+float(PCR_minus)
            SEQTot=float(SEQ_plus)+float(SEQ_minus)

            #sometimes cases of total counts being 0 occur (TODO: why?)
            if(PCRTot == 0 or SEQTot == 0):
                vcfFlagged.write("\t".join([line.rstrip(),".","."])+"\n")
                continue

            PCRFracPlus=float(PCR_plus)/PCRTot
            PCRFracMinus=float(PCR_minus)/PCRTot
            SEQFracPlus=float(SEQ_plus)/SEQTot
            SEQFracMinus=float(SEQ_minus)/SEQTot

            #simply get bias from plots
            biasSeq=int(bias_matrixSeq[current_mutation][base_before][base_after])
            biasSeqing=int(bias_matrixSeqing[current_mutation][base_before][base_after])


            if(biasSeq > 0):
                seqBiasPresent="PCR_plus=" + str(PCR_plus) + ";PCR_minus=" + str(PCR_minus)
                if(biasSeq == 1 and PCR_minus > maxNumOppositeReadsSequenceWeakBias and (PCRFracMinus > ratioVcfSequence)):
                    seqBiasPresent=0
                elif(biasSeq == 2 and PCR_minus > maxNumOppositeReadsSequenceStrongBias and (PCRFracMinus > ratioVcfSequence)):
                    seqBiasPresent=0
            if(biasSeq < 0):
                seqBiasPresent="PCR_plus=" + str(PCR_plus) + ";PCR_minus=" + str(PCR_minus)
                if(biasSeq == -1 and  PCR_plus > maxNumOppositeReadsSequenceWeakBias and (PCRFracPlus > ratioVcfSequence)):
                    seqBiasPresent=0
                elif(biasSeq == -2 and  PCR_plus > maxNumOppositeReadsSequenceStrongBias and (PCRFracPlus > ratioVcfSequence)):
                    seqBiasPresent=0
            if(biasSeqing > 0):
                seqingBiasPresent="SEQ_plus=" + str(SEQ_plus) + ";SEQ_minus=" + str(SEQ_minus)
                if(biasSeqing ==1 and SEQ_minus > maxNumOppositeReadsSequencingWeakBias and SEQFracMinus > ratioVcfSequencing):
                    seqingBiasPresent=0
                elif(biasSeqing == 2 and SEQ_minus > maxNumOppositeReadsSequencingStrongBias and SEQFracMinus > ratioVcfSequencing):
                    seqingBiasPresent=0
            if(biasSeqing < 0):
                seqingBiasPresent="SEQ_plus=" + str(SEQ_plus) + ";SEQ_minus=" + str(SEQ_minus)
                if(biasSeqing == -1 and SEQ_plus > maxNumOppositeReadsSequencingWeakBias and SEQFracPlus > ratioVcfSequencing):
                    seqingBiasPresent=0
                elif(biasSeqing == -2 and SEQ_plus > maxNumOppositeReadsSequencingStrongBias and SEQFracPlus > ratioVcfSequencing):
                    seqingBiasPresent=0
            #options to parse dir
            #if(seqBiasPresent == 0 and seqingBiasPresent == 0):
            vcfFlagged.write("\t".join([line.rstrip(),str(seqBiasPresent),str(seqingBiasPresent)])+"\n")

    return vcfFlaggedFilename



########
# MAIN #
########

# Read Parameters
parser = optparse.OptionParser()
parser.add_option('--vcfFile', action='store', type='string', dest='vcfFile', help='Specify the vcf file containing the somatic high quality SNVs.')
parser.add_option('--referenceFile', action='store', type='string', dest='referenceFile', help='Specify the filepath to the reference sequence. If this is set to "NA", then it is assumed that the VCF file contains a column with ID "SEQUENCE_CONTEXT"')
parser.add_option('--sequence_specificFile', action='store', type='string', dest='sequence_specificFile', help='Specify the filepath to sequence_specificFile. ')
parser.add_option('--sequencing_specificFile', action='store', type='string', dest='sequencing_specificFile', help='Specify the filepath to sequencing_specificFile ')

parser.add_option('--numReads', action='store', type='int', dest='numReads', help='')
parser.add_option('--numMuts', action='store', type='int', dest='numMuts', help='')
parser.add_option('--biasRatioMinimum', action='store', type='float', dest='biasRatioMinimum', help='')
parser.add_option('--biasRatioThreshold', action='store', type='float', dest='biasRatioThreshold', help='')
parser.add_option('--maxNumOppositeReadsSequenceWeakBias', action='store', type='int', dest='maxNumOppositeReadsSequenceWeakBias', help='')
parser.add_option('--maxNumOppositeReadsSequencingWeakBias', action='store', type='int', dest='maxNumOppositeReadsSequencingWeakBias', help='')
parser.add_option('--maxNumOppositeReadsSequenceStrongBias', action='store', type='int', dest='maxNumOppositeReadsSequenceStrongBias', help='')
parser.add_option('--maxNumOppositeReadsSequencingStrongBias', action='store', type='int', dest='maxNumOppositeReadsSequencingStrongBias', help='')

parser.add_option('--ratioVcf', action='store', type='float', dest='ratioVcf', help='')
parser.add_option('--biasPValThreshold', action='store', type='float', dest='biasPValThreshold', help='')

parser.add_option('--bias_matrixSeqFile', action='store', type='string', dest='bias_matrixSeqFile', help='Specify the filepath to the sequencing_specificFile ')
parser.add_option('--bias_matrixSeqingFile', action='store', type='string', dest='bias_matrixSeqingFile', help='Specify the filepath to the sequencing_specificFile ')
parser.add_option('--vcfFileFlagged', action='store', type='string', dest='vcfFileFlagged', help='Specify vcfFileFlagged.')
(options,args) = parser.parse_args()

#sys.stderr.write(' '.join(map(options,args)) + '\n')
#print(options, file=sys.stderr)

# User-specified threshold parameters default values
# sequence_specific and sequencing_specific parameters decided from plots
numReadsSequence=options.numReads
numReadsSequencing=options.numReads
numMutsSequence=options.numMuts
numMutsSequencing=options.numMuts
biasRatioMinimum=options.biasRatioMinimum
biasRatioThreshold=options.biasRatioThreshold
#biasRatioThresholdSequence=options.biasRatioThreshold
#biasRatioThresholdSequencing=options.biasRatioThreshold
# cvf thresholds
maxNumOppositeReadsSequenceWeakBias=options.maxNumOppositeReadsSequenceWeakBias
maxNumOppositeReadsSequencingWeakBias=options.maxNumOppositeReadsSequencingWeakBias
maxNumOppositeReadsSequenceStrongBias=options.maxNumOppositeReadsSequenceStrongBias
maxNumOppositeReadsSequencingStrongBias=options.maxNumOppositeReadsSequencingStrongBias
ratioVcfSequence=options.ratioVcf
ratioVcfSequencing=options.ratioVcf
biasPValThreshold=options.biasPValThreshold

# Perform Analysis
error_matrixSeq, mutation_count_matrixSeq = readSeqBiasMatrix(options.sequence_specificFile)
error_matrixSeqing, mutation_count_matrixSeqing = readSeqBiasMatrix(options.sequencing_specificFile)
bias_matrixSeq = calculateBiasMatrix(error_matrixSeq, mutation_count_matrixSeq,numReadsSequence,numMutsSequence,biasPValThreshold,biasRatioThreshold)
bias_matrixSeqing = calculateBiasMatrix(error_matrixSeqing, mutation_count_matrixSeqing,numReadsSequencing,numMutsSequencing,biasPValThreshold,biasRatioThreshold)
printErrorMatrix(bias_matrixSeq, options.bias_matrixSeqFile)
printErrorMatrix(bias_matrixSeqing, options.bias_matrixSeqingFile)

vcfFileF = markSeqBiasInVcfFile(options.vcfFile, options.referenceFile, bias_matrixSeq, bias_matrixSeqing, options.vcfFileFlagged, maxNumOppositeReadsSequenceWeakBias, maxNumOppositeReadsSequencingWeakBias, maxNumOppositeReadsSequenceStrongBias, maxNumOppositeReadsSequencingStrongBias, ratioVcfSequence, ratioVcfSequencing)
