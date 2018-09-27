#!/usr/bin/env python
#
# Copyright (c) 2018 German Cancer Research Center (DKFZ).
#
# Distributed under the MIT License (https://opensource.org/licenses/MIT).
#

# Demo for TableExtractor
# Script extracts nonsynonymous, stopgain and splicing SNVs from VCF
# Script takes VCF file produced by SNV pipeline as input

import sys
import os
import TableExtractor


fp=sys.argv[1]
outFn=sys.argv[2]
outPath=sys.argv[3]
bigger=int(sys.argv[4])-1
conf=sys.argv[4]



cNonsyn = ("EXONIC_CLASSIFICATION","in",("nonsynonymous SNV","stopgain SNV","stoploss SNV"))
cSplicing = ("ANNOVAR_FUNCTION","=","splicing")
cSomatic = ("RECLASSIFICATION","!=","lowCov_SNP_support_germline")
cGermline = ("ANNOTATION_control","!=","somatic")
cSNP = ("1K_GENOMES","1KgenomesAF_EUR",("<=",0.10))
cConf = ("CONFIDENCE",">",str(bigger))
cAnnocont = ("ANNOTATION_control","=","somatic") 
cssSomatic = [[cConf,cSomatic,cAnnocont]] # retain snvs that satisfy [condition1 and condition2] or [condition1 and condition3]
cssCoding = [[cConf,cSomatic,cNonsyn,cAnnocont],[cConf,cSomatic,cSplicing,cAnnocont]] # retain snvs that satisfy [condition1 and condition2] or [condition1 and condition3]
cssGermline = [[cGermline,cNonsyn],[cGermline,cSplicing]]



(okay,filtered) = TableExtractor.filterFile(cssSomatic,fp,os.path.join(outPath,outFn+"_somatic_snvs_conf_"+conf+"_to_10.vcf"),header='#') # the last line starting with a '#' is used to determine the header
(okay,filtered) = TableExtractor.filterFile(cssCoding,fp,os.path.join(outPath,outFn+"_somatic_coding_snvs_conf_"+conf+"_to_10.vcf"),header='#') # the last line starting with a '#' is used to determine the header
(okay,filtered) = TableExtractor.filterFile(cssGermline,fp,os.path.join(outPath,outFn+"_germline_coding_snvs_conf_"+conf+"_to_10.vcf"),header='#') # the last line starting with a '#' is used to determine the header
