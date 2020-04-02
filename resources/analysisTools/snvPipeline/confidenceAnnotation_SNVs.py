#!/usr/bin/env python
#
# Copyright (c) 2018 German Cancer Research Center (DKFZ).
#
# Distributed under the MIT License (https://opensource.org/licenses/MIT).
#

# confidenceAnnotation_SNVs.py
# Author: Jeongbin Park (based on perl code 'confidenceREAnnotation_SNVs.pl' and 'confidenceAnnotation_SNVsNoGermline.pl')
# Purpose: Annotate 'CONFIDENCE' and 'RECLASSIFICATION' to the vcf file.
# Usage: python -u confidenceAnnotation_SNVs.py --options > output.vcf
# Explanation:

import argparse
import sys
import time

from vcfparser import *  # BGZFType, LineParser, get_header_indices

def extract_info(info, keys, sep=";"):
    info_kv = {key: val for key, val in (j if len(j) == 2 else (j[0], None) for j in (i.split('=') for i in info.split(sep)))}

    if type(keys) is list:
        rtn = []
        for key in keys:
            rtn.append(info_kv.get(key, None))
            rtn = ['0' if i == 'None' else i for i in rtn]
        return rtn

    if type(keys) is str:
        rtn = info_kv.get(keys, None)
        rtn = '0' if rtn == "None" else rtn
        return rtn

def main(args):
    if args.pancanout is not None and not args.no_makehead:
        header = '##fileformat=VCFv4.1\n' \
                 '##fileDate=' + time.strftime("%Y%m%d") + '\n' \
                 '##pancancerversion=1.0\n' \
                 '##reference=<ID=' + args.refgenome[0] + ',Source=' + args.refgenome[1] + '>\n' \
                 '##center=' + args.center + '\n' \
                 '##workflowName=DKFZ_SNV_workflow\n' \
                 '##workflowVersion=1.0.0\n'
        header += '\n'.join(args.additional_header) + '\n' if len(args.additional_header) > 0 else ""
        header += '##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description="Indicates if record is a somatic mutation">\n' \
                  '##INFO=<ID=GERMLINE,Number=0,Type=Flag,Description="Indicates if record is a germline mutation">\n' \
                  '##INFO=<ID=UNCLEAR,Number=0,Type=Flag,Description="Indicates if the somatic status of a mutation is unclear">\n' \
                  '##INFO=<ID=VT,Number=1,Type=String,Description="Variant type, can be SNP, INS or DEL">\n' \
                  '##INFO=<ID=AF,Number=.,Type=Float,Description="Allele Frequency in primary data, for each ALT allele, in the same order as listed">\n' \
                  '##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP membership">\n' \
                  '##INFO=<ID=MQ,Number=1,Type=Integer,Description="RMS Mapping Quality">\n' \
                  '##INFO=<ID=1000G,Number=0,Type=Flag,Description="Indicates membership in 1000Genomes">\n' \
                  '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n' \
                  '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth at this position in the sample">\n' \
                  '##FORMAT=<ID=DP4,Number=4,Type=Integer,Description="Number of high-quality ref-forward, ref-reverse, alt-forward and alt-reverse bases">\n' \
                  '##FILTER=<ID=RE,Description="variant in UCSC_27Sept2013_RepeatMasker.bed.gz region and/or SimpleTandemRepeats_chr.bed.gz region, downloaded from UCSC genome browser and/or variant in segmental duplication region, annotated by annovar">\n' \
                  '##FILTER=<ID=BL,Description="variant in DAC-Blacklist from ENCODE or in DUKE_EXCLUDED list, both downloaded from UCSC genome browser">\n' \
                  '##FILTER=<ID=DP,Description="<= 5 reads total at position in tumor">\n' \
                  '##FILTER=<ID=SB,Description="Strand bias of reads with mutant allele = zero reads on one strand">\n' \
                  '##FILTER=<ID=TAC,Description="less than 6 reads in Tumor at position">\n' \
                  '##FILTER=<ID=dbSNP,Description="variant in dbSNP147">\n' \
                  '##FILTER=<ID=DB,Description="variant in 1000Genomes, ALL.wgs.phase1_integrated_calls.20101123.snps_chr.vcf.gz or dbSNP">\n' \
                  '##FILTER=<ID=HSDEPTH,Description="variant in HiSeqDepthTop10Pct_chr.bed.gz region, downloaded from UCSC genome browser">\n' \
                  '##FILTER=<ID=MAP,Description="variant overlaps a region from wgEncodeCrgMapabilityAlign100mer.bedGraph.gz:::--breakPointMode --aEndOffset=1 with a value below 0.5, punishment increases with a decreasing mapability">\n' \
                  '##FILTER=<ID=SBAF,Description="Strand bias of reads with mutant allele = zero reads on one strand and variant allele frequency below 0.1">\n' \
                  '##FILTER=<ID=FRQ,Description="variant allele frequency below 0.05">\n' \
                  '##FILTER=<ID=TAR,Description="Only one alternative read in Tumor at position">\n' \
                  '##FILTER=<ID=UNCLEAR,Description="Classification is unclear">\n' \
                  '##FILTER=<ID=DPHIGH,Description="Too many reads mapped in control at this region">\n' \
                  '##FILTER=<ID=DPLOWC,Description="Only 5 or less reads in control">\n' \
                  '##FILTER=<ID=1PS,Description="Only two alternative reads, one on each strand">\n' \
                  '##FILTER=<ID=ALTC,Description="Alternative reads in control">\n' \
                  '##FILTER=<ID=ALTCFR,Description="Alternative reads in control and tumor allele frequency below 0.3">\n' \
                  '##FILTER=<ID=FRC,Description="Variant allele frequency below 0.3 in germline call">\n' \
                  '##FILTER=<ID=YALT,Description="Variant on Y chromosome with low allele frequency">\n' \
                  '##FILTER=<ID=VAF,Description="Variant allele frequency in tumor < ' + str(args.newpun) + ' times allele frequency in control">\n' \
                  '##FILTER=<ID=BI,Description="Bias towards a PCR strand or sequencing strand">\n' \
                  '##SAMPLE=<ID=CONTROL,SampleName=control_' + args.pid + ',Individual=' + args.pid + ',Description="Control">\n' \
                  '##SAMPLE=<ID=TUMOR,SampleName=tumor_'+args.pid+',Individual='+args.pid+',Description="Tumor">\n'\
                  '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'
        args.pancanout.write(header)
        if not args.no_control:
            args.pancanout.write("CONTROL\t")
        args.pancanout.write("TUMOR\n")

    if args.runexome:
        args.no_punpcr = True

    idx_pcrbias = -1
    idx_seqbias = -1
    idx_pcrbias_1 = -1
    idx_seqbias_1 = -1

    for line in args.infile:
        if line[:2] == "##":
            print line.strip()
            continue

        if line[0] == "#":
            headers = list(line[1:].rstrip().split('\t'))
            fixed_headers = ["^INFO$", "MAPABILITY", "HISEQDEPTH", "SIMPLE_TANDEMREPEATS", "REPEAT_MASKER", "DUKE_EXCLUDED",
                             "DAC_BLACKLIST", "SELFCHAIN", "^CONFIDENCE$", "^RECLASSIFICATION$", "^PENALTIES$",
                             "^seqBiasPresent$", "^seqingBiasPresent$", "^seqBiasPresent_1$", "^seqingBiasPresent_1$",
                             "^seqBiasPresent_2$", "^seqingBiasPresent_2$"]
            variable_headers = { "ANNOVAR_SEGDUP_COL": "^SEGDUP$", "KGENOMES_COL": "^1K_GENOMES$", "DBSNP_COL": "^DBSNP$", }

            if args.no_control:
                variable_headers["ExAC_COL"] = "^ExAC$"
                variable_headers["EVS_COL"] = "^EVS$"
                variable_headers["GNOMAD_EXOMES_COL"] = "^GNOMAD_EXOMES$"
                variable_headers["GNOMAD_GENOMES_COL"] = "^GNOMAD_GENOMES$"
                variable_headers["LOCALCONTROL_WGS_COL"] = "^LocalControlAF_WGS$"
                variable_headers["LOCALCONTROL_WES_COL"] = "^LocalControlAF_WES$"  
            else:
                fixed_headers += [ "^INFO_control", "^ANNOTATION_control$", ]

            header_indices = get_header_indices(headers, args.configfile, fixed_headers, variable_headers)

            # create headers if they don't exist
            for optional_header in ["CONFIDENCE", "RECLASSIFICATION", ]:
                if header_indices[optional_header] == -1:
                    headers.append(optional_header)

            if args.print_info and header_indices["PENALTIES"] == -1:
                headers.append("PENALTIES")

            if args.round < 3:
                idx_pcrbias = header_indices["seqBiasPresent"]
                idx_seqbias = header_indices["seqingBiasPresent"]
            else:
                idx_pcrbias = header_indices["seqBiasPresent"] = header_indices["seqBiasPresent_2"]
                idx_seqbias = header_indices["seqingBiasPresent"] = header_indices["seqingBiasPresent_2"]
            idx_pcrbias_1 = header_indices["seqBiasPresent_1"]
            idx_seqbias_1 = header_indices["seqingBiasPresent_1"]

            if idx_pcrbias != -1 and idx_seqbias != -1:
                if args.round == 1 or args.round == 2:
                    # Update duplicating header names
                    headers[idx_pcrbias] = "seqBiasPresent_" + str(args.round)
                    headers[idx_seqbias] = "seqingBiasPresent_" + str(args.round)

                if args.round == 2 or args.round == 3:
                    if idx_pcrbias_1 == -1 or idx_seqbias_1 == -1:
                        raise ValueError("There was no column with the first bias filtering present.")
                    args.round = 2

            header_str = "#" + '\t'.join(headers)

            print(header_str)

            if args.somout is not None:
                args.somout.write(header_str + '\n')

            continue

        entries = list(line.rstrip().split('\t'))

        help = LineParser(entries, header_indices) # This will let you access values using header name

        confidence = 10 # start with maximum value and subtract something for each "bad" feature
                        # "high" = 9-10, "medium" = 6-8, low <= 5

        # in order to work with vcf files originating from mutect2 snv calling, we set confidence to 999 for these files
        # to prevent mutect2-calls from being filtered with our in house filters
        # with this means we are able to compare in house SNV calls with mutect2 SNV calls directly
        if (header_indices["CONFIDENCE"] != -1):
            confidence_string = help["CONFIDENCE"]
            if (confidence_string != ''):
                confidence_tmp = int(confidence_string)
                if confidence_tmp > 10:
                    confidence = confidence_tmp

        reasons = ""
        dbsnp_pos = None
        dbsnp_id = None

        in1KG = False
        in1KG_AF = False
        indbSNP = False

        is_commonSNP = False
        is_precious = False
        is_clinic = False
        is_repeat = False # true if SNP conicides with any of the suspicious repeat classes simple, low and satellite,
                          # which are partially redundant to other annotations
        is_STR = False # also short tandem repeats from Tandem Repeats Finder are suspicious and
                       # can conincide with other classes
        is_weird = False # coindicende with known artefact-producing regions
        if args.no_control:
            classification = "somatic" # start with default somatic
            inExAC = False
            inEVS = False
            inGnomAD_WES = False
            inGnomAD_WGS = False
            inLocalControl_WGS = False
            inLocalControl_WES = False
        else:
            # for potential re-classification (e.g. low coverage in control and in dbSNP => probably germline)
            classification = help["ANNOTATION_control"] # start with original classification

        ### For pancancer
        # genotype tumor as originally from mpileup
        gttum = entries[9].split(":")[0]
        # VAF in control will be genotype 1/1 if > 0.7 else will be 0/1
        if not args.no_control:
            gtcontr = float(extract_info(help["INFO_control"], "VAF"))
        infofield = {}
        filterfield = {}
        infofield["VT"] = "SNP" # All variants here are SNVs, no indels in this file

        # 1) external information of if these SNPs have already been found (incl. false positives from 1000 genomes!)
        # 1000 genomes
        if "MATCH=exact" in help["KGENOMES_COL"]:
            in1KG = True
            if args.no_control:
                af = extract_info(help["KGENOMES_COL"].split("&")[0], "EUR_AF")
                if af is not None and any(af > 0.01 for af in map(float, af.split(','))) > 0.01:
                    in1KG_AF = True
            infofield["1000G"] = "1000G"
        # dbSNP
        if "MATCH=exact" in help["DBSNP_COL"]:
            indbSNP = True
            infofield["DB"] = "DB"
            dbsnp_pos = entries[1]
            dbsnp_id = extract_info(help["DBSNP_COL"].split("&")[0], "ID")
            # precious!
            #INFO=<ID=CLN,Number=0,Type=Flag,Description="SNP is Clinical(LSDB,OMIM,TPA,Diagnostic)">
            #INFO=<ID=PM,Number=0,Type=Flag,Description="SNP is Precious(Clinical,Pubmed Cited)">
            if ";CLN;" in help["DBSNP_COL"]:
                is_clinic = True
            if ";PM;" in help["DBSNP_COL"]:
                is_precious = True
            if "COMMON=1" in help["DBSNP_COL"]:
                is_commonSNP = True

        if args.no_control:
            if indbSNP and is_commonSNP and not is_clinic:
                reasons += "dbSNP(NoControl)"
            if help["ExAC_COL_VALID"] and any(af > 1.0 for af in map(float, extract_info(help["ExAC_COL"], "AF").split(','))):
                inExAC = True
                infofield["ExAC"] = "ExAC"
                reasons += "ExAC(NoControl)"
            if help["EVS_COL_VALID"] and any(af > 1.0 for af in map(float, extract_info(help["EVS_COL"], "MAF").split(','))):
                inEVS = True
                infofield["EVS"] = "EVS"
                reasons += "EVS(NoControl)"

            if help["GNOMAD_EXOMES_COL_VALID"] and any(af > 0.001 for af in map(float, extract_info(help["GNOMAD_EXOMES_COL"], "AF").split(','))):
                inGnomAD_WES = True
                infofield["gnomAD_Exomes"] = "gnomAD_Exomes"
                reasons += "gnomAD_Exomes(NoControl)"
            if help["GNOMAD_GENOMES_COL_VALID"] and any(af > 0.001 for af in map(float, extract_info(help["GNOMAD_GENOMES_COL"], "AF").split(','))):
                inGnomAD_WGS = True
                infofield["gnomAD_Genomes"] = "gnomAD_Genomes"
                reasons += "gnomAD_Genomes(NoControl)"

            if help["LOCALCONTROL_WGS_COL_VALID"] and any(af > 0.01 for af in map(float, extract_info(help["LOCALCONTROL_WGS_COL"], "AF").split(','))):
                inLocalControl_WGS = True
                infofield["LocalControl_WGS"] = "LocalControl_WGS"
                reasons += "LocalControl_WGS(NoControl)"
            if help["LOCALCONTROL_WES_COL_VALID"] and any(af > 0.01 for af in map(float, extract_info(help["LOCALCONTROL_WES_COL"], "AF").split(','))):
                inLocalControl_WES = True
                infofield["LocalControl_WES"] = "LocalControl_WES"
                reasons += "LocalControl_WES(NoControl)"

        # Punish for biases round 1
        if idx_pcrbias != -1 and idx_seqbias != -1 and args.round == 1:
            if help["seqBiasPresent_VALID"] and help["seqingBiasPresent_VALID"] and not args.no_punpcr and not args.no_punseq:
                confidence -= 3
                reasons += "bias_filter_round1_PCR_and_SEQ(-3)"
                filterfield["BI"] = 1
            elif help["seqBiasPresent_VALID"] and not args.no_punpcr:
                confidence -= 3
                reasons += "bias_filter_round1_PCR(-3)"
                filterfield["BI"] = 1
            elif help["seqingBiasPresent_VALID"] and not args.no_punseq:
                confidence -= 3
                reasons += "bias_filter_round1_SEQ(-3)"
                filterfield["BI"] = 1

        # Punish for biases round 1 in round 2
        if idx_pcrbias_1 != -1 and idx_seqbias_1 != -1 and args.round == 2:
            if help["seqBiasPresent_1_VALID"] and help["seqingBiasPresent_1_VALID"] and not args.no_punpcr and not args.no_punseq:
                confidence -= 3
                reasons += "bias_filter_round1_PCR_and_SEQ(-3)"
                filterfield["BI"] = 1
            elif help["seqBiasPresent_1_VALID"] and not args.no_punpcr:
                confidence -= 3
                reasons += "bias_filter_round1_PCR(-3)"
                filterfield["BI"] = 1
            elif help["seqingBiasPresent_1_VALID"] and not args.no_punseq:
                confidence -= 3
                reasons += "bias_filter_round1_SEQ(-3)"
                filterfield["BI"] = 1

        # Punish for biases round 2
        if idx_pcrbias != -1 and idx_seqbias != -1 and args.round == 2:
            if help["seqBiasPresent_VALID"] and help["seqingBiasPresent_VALID"] and not args.no_punpcr and not args.no_punseq:
                confidence -= 3
                reasons += "bias_filter_round2_PCR_and_SEQ(-3)"
                filterfield["BI"] = 1
            elif help["seqBiasPresent_VALID"] and not args.no_punpcr:
                confidence -= 3
                reasons += "bias_filter_round2_PCR(-3)"
                filterfield["BI"] = 1
            elif help["seqingBiasPresent_VALID"] and not args.no_punseq:
                confidence -= 3
                reasons += "bias_filter_round2_SEQ(-3)"
                filterfield["BI"] = 1

        # 2) annotations of regions that cause problems: some classes of repeats from RepeatMasker track,
        # segmental duplications, (cf. Reumers et al. 2012, Nature Biotech 30:61), external blacklists, mapability
        # simple repeats and low complexity (not the same as homopolymer, but similar enough);
        # some satellites are not annotated in blacklist ...
        if any(word in help["REPEAT_MASKER"] for word in ["Simple_repeat", "Low_", "Satellite", ]):
            is_repeat = True
            confidence -= 2
            reasons += "Simple_repeat(-2)"
            filterfield["RE"] = 1
        # other repeat elements to penalize at least a bit
        elif help["REPEAT_MASKER_VALID"]:
            confidence -= 1
            reasons += "Other_repeat(-1)"
            filterfield["RE"] = 1

        # simple tandem repeats most often coincide with other bad features - do not penalize twice
        if help["SIMPLE_TANDEMREPEATS_VALID"]:
            is_STR = 1
            if not is_repeat:
                confidence -= 2
                reasons += "Tandem_repeat(-2)"
                filterfield["RE"] = 1

        # Segmental Duplications are less effective than homopolymers, short tandem repeats and microsatellites,
        # do not penality twice
        if help["ANNOVAR_SEGDUP_COL_VALID"] and not (is_repeat or is_STR):
            confidence -= 2 # bad region
            is_weird = True
            reasons += "Segmental_dup(-2)"
            filterfield["RE"] = 1

        # Duke excluded and ENCODE DAC blacklist, only consider if not already annotated as suspicious repeat
        if help["DUKE_EXCLUDED_VALID"] or help["DAC_BLACKLIST_VALID"]:
            confidence -= 3 # really bad region, usually centromeric repeats
            is_weird = True
            reasons += "Blacklist(-3)"
            filterfield["BL"] = 1

        # HiSeqDepth: regions "attracting" reads; often coincide with tandem repeats and CEN/TEL,
        # not always with low mapability
        if help["HISEQDEPTH_VALID"]:
            confidence -= 3 # really really bad region!
            is_weird = True
            reasons += "Hiseqdepth(-3)"
            filterfield["HSDEPTH"] = 1

        # Mapability is 1 for unique regions, 0.5 for regions appearing twice, 0.33... 3times, ...
        # Everything with really high number of occurences is artefacts
        # does not always correlate with the above regions
        # is overestimating badness bc. of _single_ end read simulations
        if help["MAPABILITY"] == ".":
            # in very rare cases (CEN), there is no mapability => ".", which is not numeric but interpreted as 0
            confidence -= 5
            reasons += "Not_mappable(-5)"
            filterfield["MAP"] = 1
        else:
            reduce = 0
            mapability = min(map(float, help["MAPABILITY"].split("&")))
            if mapability < 0.5:
                # 0.5 does not seem to be that bad: region appears another time in
                # the genome and we have paired end data!
                confidence -= 1
                reduce += 1

                is_weird = True # something _is_ weird already there and known SNPs might be artefacts

                if mapability < 0.4: # 3-4 times appearing region is worse but still not too bad
                    confidence -= 1
                    reduce += 1

                if mapability < 0.25: # > 4 times appearing region
                    confidence -= 1
                    reduce += 1

                if mapability < 0.1: # > 5 times is bad
                    confidence -= 2
                    reduce += 2

                if mapability < 0.05: # these regions are clearly very bad (Lego stacks)
                    confidence -= 3
                    reduce += 3

                filterfield["MAP"] = 1
                reasons += "Low_mappability(%s=>-%d)"%(help["MAPABILITY"], reduce)

        # if others have found the SNP already, it may be interesting despite low score
        # - but only if it's not a weird region.
        # if a position gets up from score 4 to 6 by giving +2 for presence in dbSNP,
        # it's very likely an artefact reported to dbSNP
        #if (in1KG or indbSNP or (args.no_control and (inExAC or inEVS or inLocalControl))) and not is_weird:
        #    confidence += 1
        #    reasons += "SNP_known_noweirdregion(+1)"

        if (args.no_control):
            # an SNV that is in dbSNP but not "clinic" or/and in 1 KG with high frequency is probably germline
            if (in1KG_AF or (indbSNP and is_commonSNP and not is_clinic) or inExAC or inEVS or inGnomAD_WES or inGnomAD_WGS or inLocalControl_WES or inLocalControl_WGS):
               classification = "SNP_support_germline"

        # 3) information from the calls and germline comparisons: coverage, strand bias, variant support, ..
        #    => can lead to reclassification
        # what we do not know: "custering", misplaced indels, in read end region (!overlapping reads)
        # (cf. Reumers et al.)
        found_dp4 = found_mq = False

        DP4tumor, MQ = extract_info(help["INFO"], ["DP4", "MQ"])
        infofield["MQ"] = "MQ=" + MQ

        # number of reference and variant reads on forward and reverse complement strand for tumor
        DP = [int(e) for e in DP4tumor.split(',')]
        # absolute coverage at position
        depthT = sum(DP) # sum in DP4 - can become zero after filter_PEoverlap.py! thus add pseudocount
        # No control: better take DP?
        trf, trr, tvf, tvr = DP

        # fraction of variant reads in tumor
        if depthT > 0:
            fr_var_tum = (tvf + tvr)/float(depthT)
        else:
            fr_var_tum = 0

        # split up in 2, to subtract 3 for tumor low coverage bc
        # that's so bad and additional penalty for low control depth
        if depthT <= 5:
            confidence -= 3
            reasons += "Tumor<5reads(-3)"
            filterfield["DP"] = 1

        if tvf < 1 or tvr < 1:	# strand bias
            confidence -= 2	# "moderate"
            if fr_var_tum <= 0.1:	# and too low variant coverage
                if args.runlowmaf:
                    confidence -= 1 # "low"
                    reasons += "Strandbias_and_MAF<=0.1_and_low_maf_filter(-3)"
                else:
                    confidence -= 3 # "low"
                    reasons += "Strandbias_and_MAF<=0.1(-5)"
                filterfield["SBAF"] = 1
            else:
                reasons+="Strandbias(-2)"
                filterfield["SB"] = 1

        if tvf + tvr <= 5:	# low variant support
            if args.runlowmaf:
                #confidence -= 0
                reasons += "ALT<=5_and_low_maf_filter(-0)"
            else:
                confidence -= 2	# "moderate"
                reasons += "ALT<=5(-2)"
            filterfield["TAC"] = 1

        if depthT >= 150 and tvf + tvr <= 2 and args.runlowmaf:
            confidence -= 3
            reasons += "Coverage_above_200_and_less_than_3_variant_reads(-3)"

        # only 1 read from each strand (might be PCR error of overlapping reads!)
        ### TODO: But overlaps are removed???
        if tvf == 1 and tvr == 1:
            if args.runlowmaf:
                confidence -= 1	# "low"
                reasons += "Only_one_ALT_per_strand_and_low_maf_filter(-1)"
            else:
                confidence -= 3	# "low"
                reasons += "Only_one_ALT_per_strand(-3)"
            filterfield["1PS"] = 1

        #if (tvf + tvr) < 6 and fr_var_tum < 0.1:  # too low variant coverage
        #    confidence -= 3	# "low"
        if tvf + tvr < 5 and fr_var_tum < 0.05: # too low variant coverage
            if args.runlowmaf:
                confidence -= 1 # "low"
                reasons += "ALT<6_and_MAF<=0.1_and_low_maf_filter(-1)"
            else:
                confidence -= 3 # "low"
                reasons += "ALT<6_and_MAF<=0.1(-3)"
            filterfield["FRQ"] = 1

        if not args.no_control:
            found_dp5 = found_tsr = False
            varsuppC = 0
            DP5control, TSR = extract_info(help["INFO_control"], ["DP5", "TSR"])

            # number of reference and variant reads on forward and reverse complement strand for control
            DP = [int(e) for e in DP5control.split(',')][:4]

            # absolute coverage at position
            depthC = sum(DP) # this does not include low base quality reads - might be 0!
            crf, crr, cvf, cvr = DP
            DP4control = ','.join(map(str, DP))
            found_dp5 = True

            # variant supporting reads in control incl. low quality
            varsuppC = int(TSR) # incl. low quality reads

            # fraction of variant reads in control: the classification called "somatic" is
            # < 1/30 of control readssupport variant
            if depthC > 0:
                fr_var_ctrl = (cvf + cvr)/float(depthC)
            else:
                fr_var_ctrl = 0

            # Punish if alternative allele frequency is low and variant is also seen in control
            # (try first with factor of tumor/control)
            # TODO: TODO: JB: Why we use the second decimal of orignal value here?
            if args.newpun != 0 and float("%.2f"%fr_var_ctrl) > 0 and classification == "somatic":
                if float("%.2f"%fr_var_tum)/float("%.2f"%fr_var_ctrl) < args.newpun:
                    confidence -= 3
                    reasons += "Alternative_allel_freq_not_" + str(args.newpun) + \
                               "_times_bigger_in_tum_than_in_control(-3)"
                    filterfield["VAF"] = 1

            # regions have an excessive number of reads in control, may not be recognized by the mapability and
            # High Seq Depth tracks (unclear if these tracks are from "normal" hg19 or 1KG reference with decoy
            # sequences!)
            # tumor could have deletion here
            if depthC > args.cutoff and not args.runexome:
                confidence -= 2
                if depthC > args.cutoff * 2: # >>300 reads at lego stack regions
                    confidence -= 2
                    reasons += "Controlcoverage>2*" + str(args.cutoff) + "(-4)"
                else:
                    reasons += "Controlcoverage>" + str(args.cutoff) + "(-2)"
                filterfield["DPHIGH"] = 1

            if depthC <= 5:
                # not _that_ bad as depthT low - will very likely become reclassified as "lowCov_SNP_support_germline"
                confidence -= 3
                reasons += "Controlcoverage<=5(-3)"
                filterfield["DPLOWC"] = 1

            if "somatic" in help["ANNOTATION_control"]:	# this includes multi_somatic
                # the overlapping reads filter may have reduced the variant read to 1 - before,
                # such SNVs were never considered
                if tvf + tvr < 2:
                    confidence -= 3	# very low!
                    reasons += "Somatic_ALT<2(-3)"
                    filterfield["TAR"] = 1

                if fr_var_ctrl > 0:	# might be germline or artefact! 1/30 reads required for calling germline
                    classification = "LQVSIG"
                    if args.controlbad:
                        confidence -= 2	# "moderate"
                        reasons += "Somatic_Control_"+str(fr_var_ctrl)+"_ALT(-2)"
                        filterfield["ALTC"] = 1
                        if fr_var_tum < 0.3:
                            confidence -= 3	# "low"
                            reasons += "Somatic_MAF<03_and_Control_"+str(fr_var_ctrl)+"_ALT(-3)"
                            filterfield["ALTCFR"] = 1

                if depthC > 0 and (varsuppC / float(depthC)) > 1/30.0:
                    # present in germline but with low qual, so likely seq error ### WAS: $varsuppC > 2.
                    # THIS WAS TOO STRINGENT FOR HIGH COVERAGE DATA
                    if classification != "LQVSIG":
                        classification = "LQVSIG"
                        if args.controlbad:
                            confidence -= 3
                            reasons += "Somatic_Control>1/30_ALT(-3)"
                            filterfield["ALTC"] = 1
                    # TODO: subtract only 2 because Rosario's exome recurrent example with LQVSIG score 7 is true! (contamination!)

                # an SNV that is in dbSNP but not "precious" or/and in 1 KG with high frequency is probably germline,
                # especially if control coverage is low (may be due to library problems)
                if in1KG or (indbSNP and not (is_precious or is_clinic)):
                    classification = "SNP_support_germline"	# missed germline?
                    if depthC < 10:
                        # Matthias' suggestion that below 10 (and not <= 5 as before)
                        # the probability to miss the variant is >95%
                        classification = "lowCov_SNP_support_germline"
                        filterfield["DB"] = 1

            # in the classification, 1/30 reads required for calling germline
            if "germline" in help["ANNOTATION_control"]:
                if fr_var_ctrl < 0.3: # too far from 50%, but that depends on coverage. I'm not convinced that qualifies for "low" / -3
                    if args.controlbad:
                        confidence -= 2 # "low"
                        reasons += "Germline_ALT<0.3(-2)"
                        filterfield["FRC"] = 1
                if in1KG or (indbSNP and not (is_precious or is_clinic)): # but this supports again - number of reads may be low!
                    classification += "_SNP_support"
                if depthC <= 10: # very probably germline due to skewed distribution at low coverage
                    classification += "_lowCov"	# => can end up as "germline_SNP_support_lowCov"

            if "LQVSIG" in help["ANNOTATION_control"]:
                # low quality variant support in germline:
                # usually a weird region; a rare case anyways
                if in1KG or (indbSNP and not (is_precious or is_clinic)):
                    classification = "SNP_support_germline"
                    filterfield["DP"] = 1
                else: # weird region
                    if args.controlbad:
                        confidence -= 2
                        reasons += "Original_LQVSIG(-2)"
                        filterfield["ALTC"] = 1

            # < 1/30 of the high quality supporting the variant AND other bases present
            # In rare cases, unclear could be LOH in tumor at a multiallelic SNP position, but we
            # cannot check that without SNV calling in control (done in LOH-BAF pipeline)
            if "unclear" in help["ANNOTATION_control"]:
                if varsuppC > 1:
                    classification = "LQVSIG"
                if in1KG or (indbSNP and not (is_precious or is_clinic)):
                    classification = "SNP_support_germline"
                    filterfield["DP"] = 1
                else: # just a weird region
                    if args.controlbad:
                        confidence -= 2
                        reasons += "Unclear_no_known_SNP(-2)"
                        filterfield["ALTC"] = 1

        if args.no_control:
            infofield["AF"] = "AF=%.2f"%fr_var_tum
        else:
            infofield["AF"] = "AF=%.2f,%.2f"%(fr_var_ctrl, fr_var_tum)

        if 'Y' in entries[0] and float(entries[5]) < 100: # CHROM and QUAL
            # SNV on chrY would have to be homozygous, too
            # het SNV on male chrX is another weird thing (unless a subpopulation had amplified chrX and mutated)
            # any SNV on chrY in female is an artefact
            confidence -= 3	# "low"
            reasons += "Y_mpileup_score<100(-3)"
            filterfield["YALT"] = 1

        #To make sure that changes in the raw filter do not influence the final result we punish them with -3
        # TODO: JB: Why we use the second decimal of orignal value here?
        if not args.runlowmaf and ((float(entries[5]) < 3 and (float("%.2f"%fr_var_tum) < 0.05 or (tvf + tvr < 5))) or (float(entries[5]) < 20 and ((tvf == 0 or tvr == 0) and (tvf + tvr <= 3)))):
            confidence -= 3
            reasons += "raw_filter_punishment(-3)"
            filterfield["FRQ"] = 1

        # to make sure that the confidence score is positive:
        if confidence < 1:
            confidence = 1

        # and not > 10
        if confidence > 10:
            confidence = 10

        if is_precious or is_clinic:
            classification += "_precious" # TODO: JB: Clinic is not acually precious. This might has to be adjusted too?

        idx_confidence = header_indices["CONFIDENCE"]
        idx_reclassification = header_indices["RECLASSIFICATION"]

        if idx_confidence == -1:
            entries.append(str(confidence))
        else:
            entries[idx_confidence] = str(confidence)

        if idx_reclassification == -1:
            entries.append(classification)
        else:
            entries[idx_reclassification] = classification

        if args.pancanout is not None:
            panout = entries[:5]
            panout.append(".")

        if dbsnp_id is not None and dbsnp_pos is not None:
            if args.pancanout is not None:
                panout[2] = dbsnp_id + "_" + dbsnp_pos
            entries[2] = dbsnp_id

        if confidence > 7:
            entries[6] = "PASS"
            if args.pancanout is not None:
                panout.append("PASS")
        else:
            if entries[6][0] == ".": # Filter column
                entries[6] = ""
            filters_line = [] if entries[6] == "" else entries[6].split(';')
            if args.pancanout is not None:
                filters_pancan = []
            for filter in ("RE","BL","DP","SB","TAC","dbSNP","DB","HSDEPTH","MAP","SBAF","FRQ","TAR","UNCLEAR","DPHIGH","DPLOWC","1PS","ALTC","ALTCFR","FRC","YALT","VAF","BI"):
                if filterfield.get(filter, 0) == 1:
                    if args.pancanout is not None:
                        filters_pancan.append(filter)
                    if not filter in filters_line:
                        filters_line.append(filter)
            if args.pancanout is not None:
                panout.append(';'.join(filters_pancan))
            entries[6] = ';'.join(filters_line)

        ### Definition of genotypes still crude, could be more precise
        if not args.no_control:
            if "germline" in help["ANNOTATION_control"] and classification != "lowCov_SNP_support_germline":
                infofield["GERMLINE"] = "GERMLINE"
                if gtcontr >= 0.75:
                    gtcontr = "1/1"
                else:
                    gtcontr = "0/1"
            elif classification == "lowCov_SNP_support_germline" or "somatic" in help["ANNOTATION_control"]:
                infofield["UNCLEAR"] = "UNCLEAR"
                if gtcontr >= 0.75:
                    gtcontr = "1/1"
                elif gtcontr >= 0.3:
                    gtcontr = "0/1"
                else:
                    gtcontr = "0/0"
            else:
                infofield["SOMATIC"] = "SOMATIC"
                gtcontr = "0/0"
                #if($gttum eq "0/0"){$gttum = "0/1";} # To change that the genotype ot the alternative is at least 0/1

        if args.pancanout is not None:
            pan_info = []
            for infooption in ("SOMATIC","GERMLINE","UNCLEAR","VT","AF","MQ","DB","1000G",):
                if infooption in infofield:
                    pan_info.append(infooption)

            panout.append(';'.join(pan_info))
            panout.append("GT:DP:DP4")

            if not args.no_control:
                panout.append(":".join( (gtcontr, str(depthC), DP4control, ) ))
            panout.append(":".join( (gttum, str(depthT), DP4tumor, ) ))
            args.pancanout.write('\t'.join(panout) + '\n')

        if args.print_info:
            idx_penalties = header_indices["PENALTIES"]
            if reasons == "":
                reasons = "."
            if idx_penalties == -1:
                entries.append(reasons)
            else:
                entries[idx_penalties] = reasons
        print('\t'.join(entries))

        if args.somout is not None and confidence > 7:
            is_out = False
            if args.no_control:
                if "somatic" in classification:
                    is_out = True
            else:
                if help['ANNOTATION_control'] == "somatic" and not "lowCov_SNP_support_germline" in classification:
                    is_out = True
            if is_out:
                args.somout.write('\t'.join(entries) + '\n')

    args.infile.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Annotate 'CONFIDENCE' and 'RECLASSIFICATION' to the vcf file.")
    parser.add_argument("-b", "--configfile", dest="configfile", nargs="?", type=argparse.FileType('r'), default=None,
                        help='Specify the config file which contains header names.')
    parser.add_argument("-w", "--nocontrol", dest="no_control", action="store_true", default=False,
                        help='Set this flag if input vcf file does not contain control information.')
    parser.add_argument("-m", "--nomakehead", dest="no_makehead", action="store_true", default=False,
                        help="Set this flag if you do NOT want to produce a pancancer conform head.")
    parser.add_argument("-o", "--pancanout", dest="pancanout", nargs="?", type=BGZFType('w'), default=None,
                        help="Outfile name including path for PanCan, will only be used if set.")
    parser.add_argument("-i", "--infile", dest="infile", nargs="?",
                        type=BGZFType('r'), default=sys.stdin,
                        help="Specify the path of vcf file. If not specified, then stdin is used instead.")
    parser.add_argument("-t", "--cutoff", type=int, dest="cutoff", default=150,
                        help="lower threshold for too high read depth in control (default 150; set to 500 " \
                             "for exomes etc.), only effective with --punishcontrol.")
    parser.add_argument("-c", "--punishcontrol", dest="controlbad", action="store_true", default=False,
                        help="lower confidence score if there are bad control events.")
    parser.add_argument("-p", "--printpenalty", dest="print_info", action="store_true", default=False,
                        help="print info on penalties into additional column 'PENALTIES' at the end.")
    parser.add_argument("-g", "--refgenome", dest="refgenome", nargs=2,
                        default=["hs37d5", "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/" \
                                           "phase2_reference_assembly_sequence/hs37d5.fa.gz", ],
                        help="reference genome used for calling ID, path (default hs37d5, ftp://ftp.1000genomes.ebi" \
                             ".ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz)")
    parser.add_argument("-z", "--center", dest="center", nargs="?", default="DKFZ",
                        help="Center (unclear if the center where the data was produced or the center of the " \
                             "calling pipeline; default DKFZ).")
    parser.add_argument("-d", "--pid", dest="pid", nargs="?", help="Patient ID (default NA).", default="NA")
    parser.add_argument("-n", "--newpunish", dest="newpun", nargs="?", type=int, default=5,
                        help="Punish if alternative is supported with low allele frequency and the alternative " \
                             "allel is present in the control (default = 5 times).")
    parser.add_argument("-s", "--addhead", dest="additional_header", nargs="+", default=[],
                        help="String with additional header line infer multiple times for multiple additional lines.")
    parser.add_argument("-f", "--somout", dest="somout", nargs="?", type=argparse.FileType('w'), default=None,
                        help="File with the high confidence somatic SNVs.")
    parser.add_argument("-a", "--round", dest="round", type=int, nargs="?", default=0,
                        help="Round of iteration for filtering (default = 0).")
    parser.add_argument("-r", "--nopunpcr", dest="no_punpcr", action="store_true", default=False,
                        help="Do not run the PCR bias filter.")
    parser.add_argument("-e", "--nopunseq", dest="no_punseq", action="store_true", default=False,
                        help="Do not run the sequencing bias filer (cannot be specified for exomes).")
    parser.add_argument("-l", "--runlowmaf", dest="runlowmaf", action="store_true", default=False,
                        help="Set this option if you want to run the low maf punishment.")
    parser.add_argument("-x", "--runexome", dest="runexome", action="store_true", default=False,
                        help="Run on exome, will turn off the high control coverage punishment " \
                             "and the PCR bias filter.")
    args = parser.parse_args()
    main(args)
