#!/usr/bin/env python
#
# Copyright (c) 2018 German Cancer Research Center (DKFZ).
#
# Distributed under the MIT License (https://opensource.org/licenses/MIT).
#

import sys
import optparse
import pysam
import numpy as np
import math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

def printErrorMatrix(error_matrix, frequency_matrix, output_filename):
	bases = ["A", "C", "G", "T"]
	error_file=open(output_filename, "w")
	possible_mutations=sorted(error_matrix.keys())
	for mutation in possible_mutations:
		error_file.write("#"+mutation+"\n")
		error_file.write("\t".join([""]+bases)+"\n")
		for base in bases:
			error_list = [str(error_matrix[mutation][base][after][0]-1)+"/"+str(error_matrix[mutation][base][after][1]-1)+";"+str(frequency_matrix[mutation][base][after]) for after in bases]
			error_file.write("\t".join([base]+error_list)+"\n")
	error_file.close()

def calculateLogSize(size, max_val, base):
	return math.log(((size/max_val)*(base-1.))+1., base)

def calculateRootedSize(size, max_val):
	if(float(size) != 0.0):
		return np.sqrt(size/max_val)
	else:
		return 0.0

def plotSquareSizeLegend(ax, colormap, min_val, max_val, max_mutation_count):
	stepsize = (max_mutation_count-1)/8.
	
	# Plot square size legend
	freq_list_legend = [[0],[0],[0],[0]]
	error_list_legend = [[0],[0],[0],[0]]
	text_list = [[""]*4, [""]*4, [""]*4, [""]*4]
	for i in range(0, 8):
		if(i==0):
			freq_list_legend[0] += [1.0]
			error_list_legend[0] += [0.0]
			text_list[0][0] = "1"
		elif(i==7):
			freq_list_legend[int(float(i)/2.)] += [float(max_mutation_count)]
			error_list_legend[int(float(i)/2.)] += [0.0]
			text_list[3][3] = str(int(max_mutation_count))
		else:
			if(i<=3):
				freq_list_legend[i] += [int(1.+i*(stepsize))]
				error_list_legend[i] += [0.0]
				text_list[i][0] = str(int(1.+i*(stepsize)))
			else:
				freq_list_legend[i-4] += [int(1.+i*(stepsize))]
				error_list_legend[i-4] += [0.0]
				text_list[i-4][3] = str(int(1.+i*(stepsize)))
	
	print "freq_list"
	print freq_list_legend
	
	hintonLegend(np.array(freq_list_legend), np.array(error_list_legend), np.array(text_list), colormap, min_val, max_val, max_weight=max_mutation_count, ax=ax)

def hinton(weight_matrix, intensity_matrix, cmap, vmin, vmax, max_weight=None, ax=None):
	"""Draw Hinton diagram for visualizing a weight matrix."""
	ax = ax if ax is not None else plt.gca()
	
	# Set colors for intensity matrix
	cm = plt.get_cmap(cmap)
	cNorm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
	scalarMap = matplotlib.cm.ScalarMappable(norm=cNorm, cmap=cm)
	intensity_colors = scalarMap.to_rgba(intensity_matrix)
	
	ax.patch.set_facecolor('gray')
	ax.set_aspect('equal', 'box')
	ax.xaxis.set_major_locator(plt.NullLocator())
	ax.yaxis.set_major_locator(plt.NullLocator())
	
	for (x,y),w in np.ndenumerate(weight_matrix):
		color = intensity_colors[x][y]
		size = 0.
		if(not(w==0)):
			size=calculateRootedSize(float(w), float(weight_matrix.max()))

		if(not(max_weight == None)):
			size = 0.
			if(not(w == 0)):
				size = calculateRootedSize(float(w), float(max_weight))

		rect = plt.Rectangle([(3-y) - size / 2, x - size / 2], size, size,
		                     facecolor=color, edgecolor=color)
		ax.add_patch(rect)
	
	plt.ylim([-1,4])
	plt.xlim(-1,4)
	ax.invert_xaxis()
	ax.invert_yaxis()

def hintonLegend(weight_matrix, intensity_matrix, text_matrix, cmap, vmin, vmax, max_weight=None, ax=None):
	"""Draw Hinton diagram for visualizing a weight matrix."""
	ax = ax if ax is not None else plt.gca()
	
	# Set colors for intensity matrix
	cm = plt.get_cmap(cmap)
	cNorm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
	scalarMap = matplotlib.cm.ScalarMappable(norm=cNorm, cmap=cm)
	intensity_colors = scalarMap.to_rgba(intensity_matrix)
	
	ax.patch.set_facecolor('gray')
	ax.set_aspect('equal', 'box')
	ax.xaxis.set_major_locator(plt.NullLocator())
	ax.yaxis.set_major_locator(plt.NullLocator())
	
	for (x,y),w in np.ndenumerate(weight_matrix):
		print str(w)
		color = intensity_colors[x][y]
		size = 0.
		if(not(w==0)):
			size = calculateRootedSize(float(w), float(weight_matrix.max()))
		
		if(not(max_weight == None)):
			size = 0.
			if(not(w==0)):
				size = calculateRootedSize(float(w), float(max_weight))
			print str(size)
		rect = plt.Rectangle([(3-y) - size / 2, x - size / 2], size, size,
					facecolor=color, edgecolor=color)
		ax.add_patch(rect)
	
	for (x,y),w in np.ndenumerate(text_matrix):
		ax.add_patch(rect)
		plt.text(3-y, x, w)


	plt.ylim([-1,4])
	plt.xlim(-1,4)
	ax.invert_xaxis()
	ax.invert_yaxis()

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

def calculateErrorMatrix(vcfFilename, referenceFilename, errorType):
	vcfFile = open(vcfFilename, "r")
	reference = None
	if(not(referenceFilename == "NA")):
		reference = pysam.Fastafile(referenceFilename)
	
	possible_mutations = ["CA", "CG", "CT", "TA", "TC", "TG"]

	# Initialize Error and Mutation Count Matrix
	error_matrix = {}
	mutation_count_matrix = {}
	for mutation in possible_mutations:
		error_matrix[mutation] = {}
		mutation_count_matrix[mutation] = {}
		possible_bases = ["A", "C", "G", "T"]
		for base_before in possible_bases:
			error_matrix[mutation][base_before] = {}
			mutation_count_matrix[mutation][base_before] = {}
			for base_after in possible_bases:
				error_matrix[mutation][base_before][base_after] = [1, 1] # Initialize with pseudo counts
				mutation_count_matrix[mutation][base_before][base_after] = 0

	header = []
	for line in vcfFile:
		if(line[0:2] == "##"):
			continue
		elif(line[0] == "#"):
			header = line[1:].rstrip().split("\t")
		else:
			split_line = line.rstrip().split("\t")

			# 23.05.2016 JB: Excluded multiallelic SNVs
			if ',' in split_line[header.index("ALT")]: continue

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

			info_list = [i.split("=") for i in split_line[header.index("INFO")].split(";")]

			# Get strand specific counts
			ACGTNacgtnPLUS = []
			ACGTNacgtnMINUS = []

			for element in info_list:
				if(element[0] == "ACGTNacgtnPLUS"):
					ACGTNacgtnPLUS = [int(i) for i in element[1].split(",")]
				elif(element[0] == "ACGTNacgtnMINUS"):
					ACGTNacgtnMINUS = [int(i) for i in element[1].split(",")]

			if (len(ACGTNacgtnPLUS)==0 or len(ACGTNacgtnMINUS)==0):
				continue
			# Count number of alternative bases
			possible_bases = ["A", "C", "G", "T", "N", "a", "c", "g", "t", "n"]
			read1_nr = ACGTNacgtnPLUS[possible_bases.index(current_mutation[1])]
			read1_r = ACGTNacgtnMINUS[possible_bases.index(current_mutation[1].lower())]
			read2_nr = ACGTNacgtnMINUS[possible_bases.index(current_mutation[1])]
			read2_r = ACGTNacgtnPLUS[possible_bases.index(current_mutation[1].lower())]

			if(errorType == "sequence_specific"):
				PCR_plus = 0
				PCR_minus = 0
				try:
					mutation_index = possible_mutations.index(current_mutation)

					PCR_plus = read1_nr + read2_r
					PCR_minus = read2_nr + read1_r

				except ValueError:
	
					current_mutation = complement(current_mutation[0])+complement(current_mutation[1])
					base_before_reverse_complement = complement(base_after)
					base_after_reverse_complement = complement(base_before)
	
					base_before = base_before_reverse_complement
					base_after = base_after_reverse_complement
	
					PCR_plus = read1_r + read2_nr
					PCR_minus = read1_nr + read2_r

				error_matrix[current_mutation][base_before][base_after][0] += PCR_plus
				error_matrix[current_mutation][base_before][base_after][1] += PCR_minus
				mutation_count_matrix[current_mutation][base_before][base_after] += 1


			elif(errorType == "sequencing_specific"):
				SEQ_plus = 0
				SEQ_minus = 0
				try:
					mutation_index = possible_mutations.index(current_mutation)

					SEQ_plus = read1_nr + read2_nr
					SEQ_minus = read1_r + read2_r

				except ValueError:
	
					current_mutation = complement(current_mutation[0])+complement(current_mutation[1])
					base_before_reverse_complement = complement(base_after)
					base_after_reverse_complement = complement(base_before)
	
					base_before = base_before_reverse_complement
					base_after = base_after_reverse_complement
	
					SEQ_plus = read1_r + read2_r
					SEQ_minus = read1_nr + read2_nr

				error_matrix[current_mutation][base_before][base_after][0] += SEQ_plus
				error_matrix[current_mutation][base_before][base_after][1] += SEQ_minus
				mutation_count_matrix[current_mutation][base_before][base_after] += 1

	return error_matrix, mutation_count_matrix

def plotErrorMatrix(error_matrix, mutation_count_matrix, output_filename, error_type, plot_title):
	possible_mutations = ["CA", "CG", "CT", "TA", "TC", "TG"]

	# Set figure properties
	figure = plt.figure(figsize=(28, 4), dpi=80)
	figure.subplots_adjust(wspace=0.1, bottom=.2, top=0.88)
	gs = gridspec.GridSpec(1, 8, height_ratios=[1], width_ratios=[1,1,1,1,1,1,1,0.2])
	colormap="PRGn"
	min_val = -3
	max_val = 3
	number_ticks_colormap = 5.

	# Calculate Maximal mutation count (necessary for normalization)
	max_mutation_count = 0
	for mutation in possible_mutations:
		for base_before in ["A", "C", "G", "T"]:
			for base_after in ["A", "C", "G", "T"]:
				if(mutation_count_matrix[mutation][base_before][base_after] > max_mutation_count):
					max_mutation_count = mutation_count_matrix[mutation][base_before][base_after]

	# Plot Square Size Legend
	ax = plt.subplot(gs[6])
	plotSquareSizeLegend(ax, colormap, min_val, max_val, max_mutation_count)
	plt.title("Number of Mutations")

	counter = 0
	for mutation in possible_mutations:
		counter += 1
		current_error_list = []
		current_frequency_list = []

		for base_before in ["T", "G", "C", "A"]:
			base_before_specific_error_list = []
			base_before_specific_frequency_list = []

			for base_after in ["A", "C", "G", "T"]:
				# Convert counts to log_2 ratios
				base_before_specific_error_list += [math.log(float(error_matrix[mutation][base_before][base_after][0])/float(error_matrix[mutation][base_before][base_after][1]),2)]
				base_before_specific_frequency_list += [mutation_count_matrix[mutation][base_before][base_after]]
			current_error_list += [base_before_specific_error_list]
			current_frequency_list += [base_before_specific_frequency_list]

		# Plotting
		ax = plt.subplot(gs[counter-1])
		hinton(np.array(current_frequency_list), np.array(current_error_list), colormap, min_val, max_val, max_weight=max_mutation_count, ax=ax)

#		ax.pcolor(np.array(current_error_list), cmap=plt.get_cmap(colormap), vmin=min_val, vmax=max_val)
		error_type_translated = plot_title

		if(mutation[0] == "C" and mutation[1] == "A"):
			plt.title(error_type_translated+"\n"+mutation[0]+"->"+mutation[1])
		else:
			plt.title(mutation[0]+"->"+mutation[1])
		if(counter == 1):
			plt.yticks([0, 1, 2, 3], ("T", "G", "C", "A"))
			plt.ylabel("preceeding")
		else:
			plt.yticks([0, 1, 2, 3], ("", "", "", ""))
		plt.xticks([0, 1, 2, 3], ["T", "G", "C", "A"])
		plt.xlabel("following")		

	ax = plt.subplot(gs[7:])
	ax.yaxis.tick_right()
	ax.imshow(np.outer(np.arange(0,1,0.01), np.ones(10)), aspect='auto', cmap=colormap, origin="lower")
	plt.xlim([0,1])
	plt.xticks([],[])
	plt.yticks(np.arange(0.*100,1.*100+((100./number_ticks_colormap)/2.),100./number_ticks_colormap), np.arange(min_val, max_val+((float(max_val)-float(min_val))/number_ticks_colormap)/2., (float(max_val)-float(min_val))/number_ticks_colormap))
	plt.title("Strand Bias")
	plt.savefig(output_filename, format="pdf")

########
# MAIN #
########

# Read Parameters
parser = optparse.OptionParser()
parser.add_option('--vcfFile', action='store', type='string', dest='vcfFile', help='Specify the vcf file containing the somatic high quality SNVs.')
parser.add_option('--referenceFile', action='store', type='string', dest='referenceFile', help='Specify the filepath to the reference sequence. If this is set to "NA", then it is assumed that the VCF file contains a column with ID "SEQUENCE_CONTEXT"')
parser.add_option('--outputFile', action='store', type='string', dest='outputFile', help='Specify the filepath to the output file (will be a png graphics file).')
parser.add_option('--errorFile', action='store', type='string', dest='errorFile', help='Specify the output file containing strand specific SNV counts.')
parser.add_option('--errorType', action='store', type='string', dest='errorType', help='Specify the error type you want to investigate (possible values are "sequencing_specific", or "sequence_specific").')
parser.add_option('--plot_title', action='store', type='string', dest='plot_title', help='Specify the title shown in the bias plot").')
(options,args) = parser.parse_args()

# Perform Analysis
error_matrix, mutation_count_matrix = calculateErrorMatrix(options.vcfFile, options.referenceFile, options.errorType)
plotErrorMatrix(error_matrix, mutation_count_matrix, options.outputFile, options.errorType, options.plot_title)
printErrorMatrix(error_matrix, mutation_count_matrix, options.errorFile)
