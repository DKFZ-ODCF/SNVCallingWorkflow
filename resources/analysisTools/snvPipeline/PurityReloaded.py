#!/usr/bin/env python
#
# Copyright (c) 2018 German Cancer Research Center (DKFZ).
#
# Distributed under the MIT License (https://opensource.org/licenses/MIT).
#

import string, sys, math, time

#Call: python PurityReloaded.py IN_VCF_FILE CONFIDENCE_COLUMN(0 based)

#Extract the chromosomal number
#def extractNum(t,prefix='chr'):
#	chrS = chrom = t[0].replace(prefix,'')
#	if chrS == 'X':
#	  return 23
#	if chrS == 'Y':
#	  return 24
#	return int(chrS)

confcol=int(sys.argv[2])
	
def prepareChromMap(prefix='chr'):
  chromMap = {}
  for i in range(1,23):
	chromMap[str(i)] = i
	chromMap["%s%i"%(prefix,i)] =i
  chromMap['X'] = 23
  chromMap['%sX'%prefix] = 23
  chromMap['Y'] = 24
  chromMap['%sY'%prefix] = 24
  return chromMap


#Compute Median Purity
def median(chrdict,outl,chromnum):
	num,med = 0.0,0.0
	for i in range(1,chromnum+1):
		if (i not in outl):
			num=num+1
			med=med+chrdict[i]
	if (num != 0):
		med = med / num
	else:
		med = 0.0
	return med


#Compute Standard deviation
def stdev(chrdict,outl,chromnum,median):
	std,num=0.0,0.0
	for i in range(1,chromnum+1):
                if (i not in outl):
                        num=num+1.0
                        std=(chrdict[i]-median)*(chrdict[i]-median)+std
	if (num > 1.0) & (std != 0):
		std=(1.0/(num-1.0))*std
	else:
		std=1.0
	return math.sqrt(std)	


def calculatSig(chromnumber):

	if (chromnumber-2 == 1):
                t = 31.82#636.6
        elif (chromnumber-2 == 2):
                t = 6-965#31.6
        elif (chromnumber-2 == 3):
                t = 4.541#12.9
        elif (chromnumber-2 == 4):
                t = 3.747#8.610
        elif (chromnumber-2 == 5):
                t = 3.365#6.869
        elif (chromnumber-2 == 6):
                t = 3.143#5.959
        elif (chromnumber-2 == 7):
                t = 2.998#5.408
        elif (chromnumber-2 == 8):
                t = 2.896#5.041
        elif (chromnumber-2 == 9):
                t = 2.821#4.781
	elif (chromnumber-2 == 10): 
		t = 2.764#4.587
	elif (chromnumber-2 == 11):
		t = 2.718#4.025
	elif (chromnumber-2 == 12):
		t = 2.681#3.930
	elif (chromnumber-2 == 13):
		t = 2.650#3.852
	elif (chromnumber-2 == 14):
		t = 2.624#3.787
	elif (chromnumber-2 == 15):
		t = 2.602#3.703
	elif (chromnumber-2 == 16):
		t = 2.583#3.686
	elif (chromnumber-2 == 17):
		t = 2.567#3.695
	elif (chromnumber-2 == 18):
		t = 2.552#3.579
	elif (chromnumber-2 == 19):
		t = 2.539#3.552
	elif (chromnumber-2 == 20):
		t = 2.528#3.527
	elif (chromnumber-2 == 21):
		t = 2.518#3.505
	elif (chromnumber-2 == 22):
		t = 2.508#3.485
	elif (chromnumber-2 == 23):
		t = 2.500#3.467
	else:
		t = 800.0

	return ((chromnumber - 1.0)/math.sqrt(chromnumber*1.0))*math.sqrt((t*t)/((chromnumber-2.0)+(t*t)))



#Performing Test
def testT(chrdict,chromnum):
	outl = []
	conti=True
	while (conti):
	        tempdict={}
	        tempdict2={}
		med=median(chrdict,outl,chromnum)
		std=stdev(chrdict,outl,chromnum,med)
	  	for i in range (1,chromnum+1):
			if (i not in outl):
				tempdict[i]=abs(chrdict[i]-med)
				tempdict2[abs(chrdict[i]-med)]=i
	  	maxi=tempdict2[max(tempdict.values())]
		
		T=tempdict[maxi]/std
		print "Testvalue: "+str(T)
	  	Z=calculatSig(chromnum-len(outl))
		print "Significants: "+str(Z)+" with a total amount of "+str(chromnum-len(outl))+" considered chr"
	  	if (T > Z):
			outl.append(maxi)
		else:
			conti=False
		
	return outl


#Performing extreme studentized deviate multiple outliner procedure
#def ESD(chrdict,chromnum):
#	outl=testT(chrdict,chromnum)
#	return outl;
	

#Check for heterozygous Mutation
#def ishet(DPfield):
#	af = allelfre(DPfield)
#	if(af<=0.6):
#		return True
#	else:
#		#print af,DPfield
#		return False

#Calculat the Allele-frequency
def allelfre(DPfield):
    try:
	zaehler=float(DPfield[2])+float(DPfield[3])
	nenner=int(DPfield[0])+int(DPfield[1])+zaehler
	if (nenner != 0):
		return (zaehler / nenner)
		
	else:
		return 1 
    except IndexError, e:
      print "\n### Error in DP field during allelfre computation:\n",DPfield
      #raise e

#Check if the Score is valid and larger than 8
def isValid(t):
#	print t
	if (t[confcol] != 'CONFIDENCE'):
		if (string.atoi(t[confcol]) > string.atoi("7")):
			return True
                else:
                        return False
        else:
                return False
            

#Extract the DP4 field of the vcf file
#def extractDP4(t):#
#	dp4=t[7].split(';DP4=')[1].split(';')[0].split(',')
#        # dp4=string.split(string.split(t[7],"DP4=")[3],",")
#        # dp4[0]=string.split(dp4[0],"=")[-1]
#        return dp4

#Extract the DP5 field of the vcf file
#def extractDP5(t):
#        dp5=string.split(string.split(t[11],";")[1],",")
#        dp5[0]=string.split(dp5[0],"=")[-1]
#        return dp5

#Calculation of GammaY
def gammay(file):	
	B=0.0
	A=0.0
	for DP4 in file:
		try:
		  #if (ishet(DP4)):
		  if allelfre(DP4) < 0.6:
		    #print "Het"
		    B+=string.atof(DP4[2])+string.atof(DP4[3])
		    A+=string.atof(DP4[0])+string.atof(DP4[1])
		except IndexError, e:
		  print e
		  print "DP4 field",DP4
	if (A+B != 0.0):
		return B/(A+B)	
	else:
		return 0.0
	
#Calculation of GammaX
def gammax(file):
	B=0.0
        A=0.0
        for DP5 in file:	
		if allelfre(DP5) < 0.6:	
		  B+=string.atof(DP5[2])+string.atof(DP5[3])
        	  A+=string.atof(DP5[0])+string.atof(DP5[1])
	if (A+B != 0.0):
                return B/(A+B)
        else: 
		return 1.0 
		
#Parse the input vcf file and calculat the "mue" values
def parseVcf(file,num):
	chrd={}
	tumor={}
	ref={}
	for i in range (1,num+1):
	  tumor[i] = []
	  ref[i] = []
	chromMap = prepareChromMap(prefix='chr')  # dicionary to map chromosome symbols to integers
	
	s=time.time()
	print "Start processing"
	if len(file) > 3 and file[-3:] == ".gz":
		import gzip
		vcf = gzip.open(file,"r")  
		l=vcf.readline()
	else:
		vcf = open(file,"r")
		l=vcf.readline()
	
	while (l!= ""):
		t=l.split('\t')
		if (t[0][0] != "#") and isValid(t):
			i = chromMap[t[0]]
			if (t[12]=="germline"):
			  #DP5=string.split(string.split(t[11],";")[1],",")
			  #ref[i].append(string.split(DP5[0],"=")[-1])
			  ref[i].append(t[11].split(';DP5=')[1].split(';')[0].split(','))
			if (t[12]=="somatic"):
			  tumor[i].append(t[7].split(';DP4=')[1].split(';')[0].split(','))
		l=vcf.readline()
	vcf.close()
	
	for i in range (1,num+1):
		refLs = ref[i]
		tumorLs= tumor[i]
		print "Germline mutations: "+str(len(refLs))
		print "Somatic mutations: "+str(len(tumorLs))
		gy=gammay(tumorLs)
		gx=gammax(refLs)
		print "Gammax: "+str(gx)
		print "Gammay: "+str(gy)
		chrd[i]=gy/gx
	print"Requested time: "+str(time.time()-s)  
	return chrd

	

def main():
	#Reading parameters
	inputfilename= sys.argv[1]
	amountofchr= 24#string.atoi(sys.argv[2])
#	confcol=string.atoi(sys.argv[2])
	chrpurdict= {}
	outlierl= []
	#Open files and calculate gammas and chromosomal purity
	#Do work for all "normal" chromosomes
	start=time.time()
	chrpurdict = parseVcf(inputfilename,amountofchr)

	#Output of chromosomal purity
	print "Estimated purity for every single chromosome:"
	for i in chrpurdict.items():
		print i

	print "Median purity before ESD:"
	m=(median(chrpurdict,outlierl,amountofchr))
	print m
	print "Standard deviation before ESD:"
	print str(stdev(chrpurdict,outlierl,amountofchr,m))
	#Apply ESD
	outlierl=testT(chrpurdict,amountofchr)
	print "Outlier:"+str(outlierl)
	print "Median purity after ESD"
 	print str(median(chrpurdict,outlierl,amountofchr))
	print "Standard deviation after ESD"
	print str(stdev(chrpurdict,outlierl,amountofchr,median(chrpurdict,outlierl,amountofchr)))
	
	#Output of our results
	print "Complet time: "+str(time.time()-start)
main()
