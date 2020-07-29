# author: Oliver Alka
# date: 29.05.2019
#
# Tool to extract decoy information from rerooted trees (passatuto) and parse them into the original
# assay library structure from OpenMS::AssayGeneratorMetabo
#

# packages 
import os
import re
import numpy
import pandas as pd
pd.options.mode.chained_assignment = None
import tempfile
from pyopenms import *

# constant
mass_ch2 = 14.0156500641

# functions
def addCH2(x):
	return(x + mass_ch2)

# create validated temporary file (.tsv) from OpenMS assay library input (.tsv, .pqp, .traml)
def fillTmpTSVWithValidTargetedExp(openmslib, tmpfile):
	# check file extension, validate, convert to tsv (if necessary)
	filename, extension = os.path.splitext(openmslib)
	targeted_exp = TargetedExperiment()

	if extension == '.pqp':
		TransitionPQPFile().convertPQPToTargetedExperiment(openmslib.encode(), targeted_exp, False)
	elif extension == ".traML" or extension == ".TraML" or extension == ".traml":
		TraMLFile().load(openmslib.encode(), targeted_exp)
	else:
		filetype = FileTypes().nameToType('TSV')
		TransitionTSVFile().convertTSVToTargetedExperiment(openmslib.encode(), filetype, targeted_exp)

	# check validity of OpenMS::TargetedExperiment
	TransitionTSVFile().validateTargetedExperiment(targeted_exp)
	TransitionTSVFile().convertTargetedExperimentToTSV(tmpfile.name.encode(), targeted_exp)

	return tmpfile

# extract informations from decoy trees  
def extractInfoFromRerootedTrees(treeextracts, treepath, min_fragment_mz, max_fragment_mz):
	# initialize class
	transition = PossibleTransition()
	# extraction patterns for metabolite information: 1Cymiazolhydrochloride1377_C12H14N2S_M+H+.txt
	basename = os.path.basename(treepath)
	patternsf = '_(\w.*)_' # C12H14N2S
	patternad = '_(M.*)\.' # M+H+
	patternna = '\d(\w.*\d)\d_' # 1Cymiazolhydrochloride1377
	patternall = [patternsf, patternad, patternna]

	# extract information filename (sumformula, adduct, compoundname)
	for pattern in patternall:
		extract = re.search(pattern, basename)
		if str(pattern) == str(patternsf) and extract:
			formula = extract.group(1)
			transition.sumformula = formula
		elif str(pattern) == str(patternad) and extract:
			adduct = extract.group(1)
			transition.adduct = adduct
		elif str(pattern) == str(patternna) and extract:
			cname = extract.group(1)
			transition.compoundname = ''.join([i for i in cname if not i.isdigit()])

	# parse metabolite file and extract precursormass, mass, intensity
	tree = open(treepath, "r")
	lines = tree.readlines()
	i = 0
	maxint = 0.0
	ms2sumint = 0.0
	mass=[]
	intensity=[]
	while i < len(lines):
		line = lines[i]
		if "CH$EXACT_MASS" in line:
			pmass = lines[i].split(" ")[1].strip()
			transition.precursormass = pmass
		if "PK$PEAK" in line:
			i+=1
			while "//" not in lines[i]:
				line = lines[i]
				peak_list=line.split(" ")
				mz=float(peak_list[2])
				mz=round(mz,3)
				absint=round(float(peak_list[3].split("\n")[0]),3)
				mass.append(mz)
				intensity.append(absint)
				ms2sumint += absint 
				if absint > maxint:
					maxint = absint
				i+=1
			# calculate relative intensity as in original assay library.
			if maxint != 0:
				intensity[:] = [x / maxint for x in intensity]
		else:
			i+=1

	# Remove mass and intensity of peaks below the min_fragment_mz
	i = 0
	index_below_min_threshold = []
	for i in range(0,len(mass)):
		if mass[i] < min_fragment_mz:
			index_below_min_threshold.append(i)
	for index in sorted(index_below_min_threshold, reverse=True): # remove the highest index first
		del mass[index]
		del intensity[index]
		
	# Remove mass and intensity of peaks higher the max_fragment_mz
	i = 0
	index_higher_max_threshold = []
	for i in range(0,len(mass)):
		if mass[i] > max_fragment_mz:
			index_higher_max_threshold.append(i)
	for index in sorted(index_higher_max_threshold, reverse=True): # remove the highest index first
		del mass[index]
		del intensity[index]

	# sort by intensity -> uses highest intensity peaks first 
	mass = numpy.array(mass)
	intensity = numpy.array(intensity)
	ints = intensity.argsort()[::-1][:len(intensity)]
	sort_mass = mass[ints]
	sort_int=intensity[ints]

	transition.metaboliteid = transition.compoundname + transition.adduct
	transition.mass = sort_mass
	transition.intensity = sort_int
	transition.ms2sumint = ms2sumint
	treeextracts.append(transition)
	tree.close()

	return treeextracts

# make unique list based on the highest sum intensity MS2 
# of the same compound and adduct combination
def makeUniqueCompoundAdductHighestMS2Int(treeextracts):
	highest_ms2int = {}
	# find maximum intensity for a specific entry (identifier)
	for entry in treeextracts:
		current_id = entry.metaboliteid 
		current_int = entry.ms2sumint
		highest_ms2int[current_id] = current_int
		for metabolite in treeextracts:
			if current_id == metabolite.metaboliteid and current_int < metabolite.ms2sumint: 
				highest_ms2int[metabolite.metaboliteid] = metabolite.ms2sumint
				current_id = metabolite.metaboliteid
				current_int = metabolite.ms2sumint

	# unique tree extracts 
	unique_extracts = []
	for key,value in highest_ms2int.items():
		for entry in treeextracts:
			if key == entry.metaboliteid and value == entry.ms2sumint:
				unique_extracts.append(entry)

	return unique_extracts
	
# class
class PossibleTransition:
	metaboliteid = "" 
	compoundname = ""
	sumformula = ""
	adduct = ""
	decoy = 1
	precursormass = 0.0
	mass = []
	intensity = []
	ms2sumint = 0.0

def main(openmsassaylib, rerootedtrees, decoyassaylib, addms1shift, min_fragment_mz, max_fragment_mz):
	print('Create temporary file')
	# create validated temporary file 
	tmpfile = tempfile.NamedTemporaryFile(mode='w+t', suffix='.tmp', delete=False)
	fillTmpTSVWithValidTargetedExp(openmsassaylib, tmpfile)

	# check output file extension
	base, extension = os.path.splitext(decoyassaylib)
	if extension == ".tsv" or extension == ".pqp":
		output_type = extension
	else:
		print('Output type is invalid! Please specify a valid output_type (.tsv, .pqp).')
		print('If output_type is not defined .pqp will be used automatically')
		output_type = ".pqp"
	
	output = base + output_type
	
	print('Extract information from rerooted trees')

	treeextracts = []
	for entry in os.listdir(rerootedtrees):
		if entry.endswith(".txt"):
			treepath = rerootedtrees + entry
			extractInfoFromRerootedTrees(treeextracts, treepath, min_fragment_mz, max_fragment_mz)

	# unique list based on the highest sum intensity MS2 
	# of the same compound and adduct combination
	
	unique_extracts = makeUniqueCompoundAdductHighestMS2Int(treeextracts)

	print('Read assay library input')
	
	# read input assay library and get unique TransitionGroupID
	# this step is relevant to split the dataframe into smaller instances,
	# which is then filled step by step depending on the existing
	# possibleTransition class (possible decoys)
	
	assaylib = pd.read_csv(tmpfile, sep = '\t')
	unique_tgid = assaylib.TransitionGroupId.unique()

	# split dataframe by TransitionGroupID 
	d_lib = dict()
	d_lib_decoy = dict()
	for element in unique_tgid:
		d_lib[element] = assaylib[assaylib['TransitionGroupId'] == element]
		d_lib_decoy[element] = assaylib[assaylib['TransitionGroupId'] == element]
		
	print('Generate decoys')
	
	# In general the rerooted fragmentation tree is used for decoy generation (MS2). 
	# For MS1, since it is used for scoring from OpenSWATH, the initial mass of the 
	# precursor has to be changed. Since for most organic compounds an addition of
	# CH2 is feasible, we added the CH2 mass to the precursor. 
	
	# Additionally, if either know rerooted tree was generated, or the generate fragments 
	# are similar to the one from the target, CH2 was also added to the fragment masses. 
	# This only happens in around 5% of the cases. 
	
	# Iterate over entries in list(PossibleTransitions)
	# find overlap of CompoundName, SumFormula and Adduct
	# add decoy, rows, based on n highest intensity peaks in spectrum 
	# n highest depends on the number of transitions per metabolite in the original
	# assay library
	
	# if tree was rerooted with the corresponding adduct
	count_reroot_adduct = 0
	count_reroot_adduct_add_individual_CH2 = 0 
	for key, value in d_lib_decoy.items():
		n_transitions = len(value["CompoundName"])
		db_name = value['CompoundName'].unique().tolist()[0]
		db_adduct = value["Adducts"].unique()
		for entry in unique_extracts:
			if str(entry.compoundname) == db_name and str(entry.adduct) == db_adduct and db_name.find("_decoy") == -1:
				value['TransitionGroupId'] = value['TransitionGroupId'].astype(str) + '_decoy'
				value['TransitionId'] = value['TransitionId'].astype(str) + '_decoy'
				value['CompoundName'] = value['CompoundName'].astype(str) + '_decoy'
				value['Decoy'] =  1
				value['Annotation'] = 'NA'
				# issue: values are overwritten with decoy information
				# if there are not enough decoy intensities/fragments
				# replace the missing ones with 0.0 - do the same for ProductMz
				if (len(entry.intensity) < n_transitions):
					n_missing = n_transitions - (len(entry.intensity))
					additional_transitions = [0]*n_missing
					value['LibraryIntensity'] = entry.intensity[:n_transitions].tolist() + additional_transitions
				else:
					value['LibraryIntensity'] = entry.intensity[:n_transitions]
				if addms1shift:
					value['PrecursorMz'] = addCH2(value['PrecursorMz'])
				
				# if the fragment masses are too close to each other (target and decoy) - remove them
				# and add CH2 to the fragments.
				if bool(set(list(map(int,value['ProductMz'].tolist()))).intersection(list(map(int,entry.mass[:n_transitions])))):
					# careful currently working with int values instead of actual ones 
					overlapping_values = list(set(list(map(int,value['ProductMz'].tolist()))).intersection(list(map(int,entry.mass[:n_transitions]))))
					overlapping_index=[]
					for o_value in overlapping_values:
						overlapping_index.append(list(map(int,entry.mass[:n_transitions])).index(o_value))
					current_decoy_mass = list(entry.mass[:n_transitions])
					# replace by adding ch2 mass to overlapping fragment
					for o_index in overlapping_index:
						current_decoy_mass[o_index] = addCH2(current_decoy_mass[o_index]) 
					value['ProductMz'] = current_decoy_mass
					count_reroot_adduct_add_individual_CH2 +=1
					continue;
				elif (len(entry.mass) < n_transitions):
					n_missing = n_transitions - (len(entry.mass))
					additional_transitions = [0]*n_missing
					value['ProductMz'] = entry.mass[:n_transitions].tolist() + additional_transitions
				else:
					value['ProductMz'] = entry.mass[:n_transitions]
				count_reroot_adduct += 1

	# if the method above can not be applied, but a full list of decoys should be generated
	# the remaining ones are a shifted by a -CH2 mass (Precursor and Fragments) 
	count_shifted = 0
	for key, value in d_lib_decoy.items():
		db_name = value['CompoundName'].unique().tolist()[0]
		# if "_decoy" is not found in the name proceed:
		if db_name.find("_decoy") == -1:
			value['TransitionGroupId'] = value['TransitionGroupId'].astype(str) + '_decoy'
			value['TransitionId'] = value['TransitionId'].astype(str) + '_decoy'
			value['CompoundName'] = value['CompoundName'].astype(str) + '_decoy'
			value['Decoy'] =  1
			value['Annotation'] = 'NA'
			current_fmz = value['ProductMz'].tolist()
			fmz = [addCH2(x) for x in current_fmz]
			value['ProductMz'] = fmz
			if addms1shift:
				value['PrecursorMz'] = addCH2(value['PrecursorMz'])
			count_shifted += 1

	# create df from dicts
	print("Create dataframe") 
	new_df = pd.DataFrame()
	for key, value in d_lib.items():
		new_df = new_df.append(value)
	for key, value in d_lib_decoy.items():
		new_df = new_df.append(value)

	new_df.to_csv(tmpfile.name, sep = '\t') 

	# convert tsv to TargetedExperiment - to pqp
	print("Export data")

	targetedexp = TargetedExperiment()
	transitiontsv = TransitionTSVFile()

	# validate decoy assay library 
	filetype = FileTypes().nameToType('TSV')
	transitiontsv.convertTSVToTargetedExperiment(tmpfile.name.encode(), filetype, targetedexp)
	transitiontsv.validateTargetedExperiment(targetedexp)

	if output_type == '.tsv':
		transitiontsv.convertTargetedExperimentToTSV(output.encode(), targetedexp)
	else:
		transitionpqp = TransitionPQPFile()
		transitionpqp.convertTargetedExperimentToPQP(output.encode(), targetedexp) 

	# remove temporary file 
	tmpfile.close()
	os.remove(tmpfile.name)

	print("# of decoys generated by rerooted trees (matching adduct): " + str(count_reroot_adduct))
	print("# of decoys generated by rerooted trees with one or mutiple overlapping fragments enhanced by adding CH2: " + str(count_reroot_adduct_add_individual_CH2))
	print("# of decoys generated by CH2 mass shift: " + str(count_shifted))
	print("Done")
