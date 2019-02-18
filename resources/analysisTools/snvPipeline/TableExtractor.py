#
# Copyright (c) 2018 German Cancer Research Center (DKFZ).
#
# Distributed under the MIT License (https://opensource.org/licenses/MIT).
#

# Author: Lars Feuerbach
# Extractor for tab seperated files

import sys, string

def defineExampleFilter():
  """Returns a list of tuples that represent the filter conditions."""
  ls  = []
  ls.append( (0,"=","chr1") ) # first collum has entry chr1
  ls.append( ("SEGDUP",'=','.') ) # column with header "SEGDUP" has value '.'
  ls.append( ("CONFIDENCE",'>',6) ) # column with header "CONFIDENCE" has a value greater than 6
  
  #c1 = (0,"=","DP=10")
  #ls.append( ("INFO","split",(';',c1))) # split the text in the column INFO by ';' and apply condition c1 to it
  
  ls.append(("1K_GENOMES","1KgenomesAF",(">=",0.10)))
  ls.append(("1K_GENOMES","1KgenomesAF",("<=",0.30)))
  return ls


def getHeader(f,header,sep,fileend):
  """Retrun the last line of the header"""
  if type(header) == type(1):
    for i in range(0,header-1):
      f.readline()
    return string.split(f.readline()[:-len(fileend)],sep) # return list of column names
	
  else: # if header is a character
	
	oldline = f.readline()
	if oldline[0] != header:
	  f.seek(0,0) # jump to start of file
	  return False
	pos = f.tell() # absolute position of the files read pointer
	line = f.readline()
	while line:
	  #print line
	  #print pos
	  if line[0] != header: # the last line was the header
	    f.seek(pos,0) # jump to first position after the header
	    return string.split(oldline[:-len(fileend)],sep) # return list of column names
	  oldline = line
	  line = f.readline()
	  pos = f.tell() # absolute position of the files read pointer
	# the file is empty or all lines start with header character
	return False




def parseConditions(conditions,header):
  """Checks the format of the <conditions> and the consitency with the file header """
  
  try:
    if type(conditions[0]) != type([1]): # if conditions is not a list of lists, reformt the variable
      conditions = [conditions]
  except Exception, e:
    print e, conditions , "\n Conditions are empty or do not met required format."
    sys.exit(1)
  
  for cs in conditions:
    # Check if column ids can be resloved to column numbers where necessary 
    for i in range(0,len(cs)):
      (colId,operation,parameter) = cs[i]
      if type(colId) == type(''):
	if header == False:
	  print "Error in condition %s\n"%str(cs[i]), "File has no header, thus Columnname %s can not be attributed to a column number"%colId
	  sys.exit(1)
	try:
	  cs[i] = (header.index(colId),operation,parameter)
	except ValueError:
	  print "Error in condition %s\n"%str(cs[i]), "Columnname %s is not part of file header:\n%s"%(colId,header)
	  sys.exit(1)
      
      if type(colId) == type([]):
	for j in range(0,len(colId)):
	  if type(colId[j]) == type(''):
	    if header == False:
	      print "Error in condition %s\n"%str(cs[i]), "File has no header, thus Columnname %s can not be attributed to a column number"%colId
	      sys.exit(1)
	    try:
	      colId[j] = header.index(colId[j])
	    except ValueError:
	      print "Error in condition %s\n"%str(cs[i]), "Columnname %s is not part of file header:\n%s"%(colId,header)
	      sys.exit(1)
  
  return conditions

def conditionOkay(c,cols):
  """Checks if condition <c> is satisfied by <cols> and <header> """
  (colId,operation,parameter) = c
  
  # extract values that are relevant for filter operation
  value = None
  try:
    if type(colId) in (type(1),type('')): # single colId
      value = cols[colId]
    
    elif type(colId) in (type((1,1)),type([])): # multiple colIds
      value = []
      for cId in colId:
	value.append(cols[cId])
  
  except IndexError, e:
    print cols
    print c
    print colId, "is higher then row length ", len(cols)
    raise e
  
  # execute operations
  
  # atomic operations
  if operation == "=": return value == parameter
  if operation == "!=": return value != parameter
  try:
    if operation == "==": return float(value) == float(parameter)
    if operation == ">": return float(value) > float(parameter)
    if operation == "<": return float(value) < float(parameter)
    if operation == "<>": return float(value) <> float(parameter)
    if operation in [">=","=>"]: return float(value) >= float(parameter)
    if operation in ["<=","=<"]: return float(value) <= float(parameter)
  except ValueError, vErr:
    if value != '.' and value != '-NA-': 
      print value, parameter
      raise vErr
  if operation == "in": return (value in parameter)
  
  
  # recursive operations
  if operation == "ins": # if one list item is in value
    for item in parameter:
	#print item,"in",value
	if item in value: return True
    return False
  if operation == "split":
    newCols = string.split(value,parameter[0])
    return conditionOkay(parameter[1],newCols)
  
  # special operations
  # execute a condition on the allele frquency from 1000genomes
  # Example: ("1K_GENOMES","1KgenomesAF",("<=",0.10))
  if operation == "1KgenomesAF": 
    af = "0.0"
    if value != '.' :
      i = value.index("AF=")
      af = value[i+3:i+7]  
    return conditionOkay((0,parameter[0],parameter[1]),[af]) # execute operation specified by <parameter[0]> with parameter <parameter[1]> on allele frequence
  
  if operation == "1KgenomesAF_EUR": 
    af = "0.0"
    if value != '.' :
      i = 0
      try:
	i = value.index("EUR_AF=")
	af = value[i+7:i+11]
      except ValueError:
	i = value.index("AF=")
	af = value[i+3:i+7]  
    return conditionOkay((0,parameter[0],parameter[1]),[af]) # execute operation specified by <parameter[0]> with parameter <parameter[1]> on allele frequence
  
  # Returns True if SNV is embeded in sequence context <parameter>
  # (("REF","SEQUENCE_CONTEXT"),'seqcontext',"CG") retuns all SNVs that disrupt a CG
  # (("VAR","SEQUENCE_CONTEXT"),'seqcontext',"CG") retuns all SNVs that generate a CG
  if operation == "seqcontext":
      ref = value[0]
      context = string.split(value[1],',') # seperate string "NNNNNNNNN,NNNNNNNNNN" in leading and tailing context
      patternLength = len(parameter)
      #print context
      #print parameter
      #print context[0][-(patternLength-1):]
      
      return context[0][-(patternLength-1):]+ref == parameter or ref+context[1][:(patternLength-1)] == parameter # ACGT,ACGT and C and CG -> if TC == CG or CA == CG: return True

def conditionsOkay(cs,cols):
  """If one of the conditions in <cs> is violated according to the value list <cols> the function returns False (AND linked conditions) """
  for c in cs:
    try:
      if not conditionOkay(c,cols): return False
    except IndexError, e:
      print "!!! Line is shorter than expected, maybe input file is corrupted !!!"
      return False
      #sys.exit(1)
  return True


def countFile(conditions,names,fp,header='#',sep='\t',lineend='\n',debug=False):
  """
  Count number of lines that meet <conditions>. Each counts for each conditions are stored under corresponding name <names>.
  """
  f = open(fp,'r')
  #if op: out = open(op,'w')
  if len(conditions) != len(names):
    print "Problem in count File: conditions and names have to be of equal length"
    sys.exit(1)
  #out = sys.stdout
  cnt = {}
  for name in names:
    cnt[name] = 0
  
  # retrive header and set file pointer to first line after header
  print "@Header", header
  if header: #and len(header) < 2:
    header = getHeader(f,header,sep,lineend)
    #if header and op: out.write(sep.join(header)+lineend)
  
  if debug: print "Header",header
  # check format and consitency of conditions with header
  if debug: print "\nInitial conditions:\n", conditions
  conditions = parseConditions(conditions,header)
  if debug: print "\nResolved conditions:\n", conditions
  
  for line in f:
    cols = string.split(line[:-len(lineend)],sep) # generate a list that contain all values of this row seperated by columns
    for i in range(0,len(conditions)):
      if conditionsOkay(conditions[i],cols) :
		cnt[names[i]] += 1
  #f.close()
  
  s = ""
  for name in names:
    s += "\t%i"%cnt[name]
  return s

def countVector(conditions,name,fp,header='#',sep='\t',lineend='\n',debug=False):
  """
  Count number of lines that meet <conditions>. Each counts for each conditions are stored under corresponding name <names>.
  """
  f = open(fp,'r')
  cnt = {}
  
  print "@File",fp
  # retrive header and set file pointer to first line after header
  print "@Header", header
  if header: #and len(header) < 2:
    header = getHeader(f,header,sep,lineend)
    #if header and op: out.write(sep.join(header)+lineend)
  nameCol = 0
  #infoCols = []
  try:
	nameCol = header.index(name)
  except IndexError:
	  print name, "is not contained in header:",header
  if debug: print "Header",header
  # check format and consitency of conditions with header
  if debug: print "\nInitial conditions:\n", conditions
  conditions = parseConditions(conditions,header)
  if debug: print "\nResolved conditions:\n", conditions
  
  for line in f:
    cols = string.split(line[:-len(lineend)],sep) # generate a list that contain all values of this row seperated by columns
    for i in range(0,len(conditions)):
      if conditionsOkay(conditions[i],cols) :
			cnt[cols[nameCol]] = 1
      else:
			cnt[cols[nameCol]] = 0
  f.close()
  
  return cnt

def filterFile(conditions,fp,op,header='#',sep='\t',lineend='\n',debug=True):
  """Filter a text file table according to given list of <conditions>.
     A list of condition lists. Conditions in the inner lists are compared by logical ANDs while the other conditions are evaluated by logical ORs.
     <fp> is the filepath of the input file.
     <op> the filepath of the output filepath
     <header> is either False, or the symbol that marks the header line or the numbe of the line that contains the header.
     <sep> is the character that is used to seperate the columns of the table.
     <lineend> is the character string that marks the end of each line.
  """
  f = None
  if fp=="-":
    f = sys.stdin
  elif len(fp) > 3 and fp[-3:] == ".gz":
	import gzip
	f = gzip.open(fp,'r')  
  else:
    f = open(fp,'r')
  if op: out = open(op,'w')
  #out = sys.stdout
  cntOkay = 0
  cntFiltered = 0
  # retrive header and set file pointer to first line after header
  print "@Header", header
  if header: # and len(header) < 2:
    header = getHeader(f,header,sep,lineend)
    if header and op: out.write(sep.join(header)+lineend)
  if debug: print "Header",header
  # check format and consitency of conditions with header
  if debug: print "\nInitial conditions:\n", conditions
  conditions = parseConditions(conditions,header)
  if debug: print "\nResolved conditions:\n", conditions
  oldLine = None
  for line in f:
    
    if line == oldLine: continue # skip lines that are contained multiple times
    
    cols = string.split(line[:-len(lineend)],sep) # generate a list that contain all values of this row seperated by columns
    for cs in conditions:
      if conditionsOkay(cs,cols) :
	if op: out.write(line)
	cntOkay += 1
	break # only one of the cs in conditions has to be met (OR linked)
    cntFiltered += 1
    oldLine = line
  if op: out.close()
  #f.close()
  print "\n%s -> %s\n %i okay\t%i filtered\t%i total\t%.4f %%\n"%(fp,op,cntOkay,cntFiltered,cntOkay+cntFiltered,(100.0*cntOkay)/(cntOkay+cntFiltered))
  return (cntOkay,cntFiltered)

## Test call
##filterFile(defineExampleFilter(),"/home/feuerbac/LSDF_Prostate/results_per_pid/1111492/mpileup/snvs_1111492.vcf","")
#conditions = [[("score",">",0.1)]]
#fp = "/home/feuerbac/Data/ICGC_Prostate/Integration/SmallVsLarge/Scored/ICGC_PCA018_Integrated_2013124.tab"
#cnt = countVector(conditions,"gene_name",fp,debug = 5)
##cnt = countVector([[("score",">",0.1)]],"gene_name","/home/feuerbac/Data/ICGC_Prostate/Integration/SmallVsLarge/Scored/ICGC_PCA018_Integrated_2013124.tab",header='#',sep='\t',lineend='\n',debug=False)
##for i in range(0,100):
##	key = cnt.keys()[i]
##	print key,cnt[key]
#print cnt["JUN"]
