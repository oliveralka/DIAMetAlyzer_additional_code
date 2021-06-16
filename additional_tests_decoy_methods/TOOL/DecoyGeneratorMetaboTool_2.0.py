# DecoyGeneratorMetabo
#
# Tool to generate decoys for an metabolomics assay library using different methods
#
# author: Oliver Alka
# date: 05.0.2019 
#
# methods: 
# 1: rtperm - retention time permutation (same precursor and fragments - other RT) 
# 2: mz (same precursor and fragments - other mz or other window))
# 3: linmassperm - same precursor mz, but different fragments (linear mass shift)
# 4: swathwindow - precursor_mz to other mass window
# 5: fragdb - use the same precursor, but fragments with lower mz as prec from
#             all the fragment masses existing in the assay annotatet library.

# The passatuto fragmentation tree based method can be used using the DIAMetAlyzer 
# KNIME workflow 


# packages
import os
import click
import random
from pyopenms import *
import pandas as pd
pd.options.mode.chained_assignment = None
import tempfile
from pyopenms import *
from collections import Counter

# constant
mass_ch2 = 14.015650

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

# click command
@click.command()
@click.option('--assay_library', '-as', envvar = 'assay_library', multiple = False, type = click.Path(), help = 'Assay library (AssayGeneratorMetabo - .tsv only) ')
@click.option('--swath_windows', '-sw', envvar = 'swath_windows', multiple = False, type = click.Path(), help = 'Swath windows - .txt')
@click.option('--output', '-out', envvar = 'output', multiple = False, type = click.Path(), help = 'Path to save the output file')
@click.option('--method', '-m', envvar = 'method', multiple = False, default="linmzperm", show_default=True, type = click.Choice(["linmzperm","rtperm","shufflefrag","sw_perm","fragdb"]), help = 'Methods for decoy creation')
@click.option('--rt_mindist', '-rtmd', envvar = 'rt_mindist', multiple = False, default=80.0, show_default=True, type = float, help = 'Minimal distance from current retention time - has to be set higher than rt distance used in OpenSWATH')
@click.option('--retentiontime', '-rt', envvar = 'retentiontime', multiple = False, default=950, show_default=True, type = float, help = 'Retention time of gradient used e.g. 950 seconds (retention time between 0 and 950 seconds)')
@click.option('--addCH2ToPrec', '-ch2prec', envvar = 'addch2toprec', is_flag=True, help="Shifts precursor mass by CH2 (works methods on MS2 level; linmzperm, shufflefrag, fragdb)")

def main(assay_library, swath_windows, output, method, rt_mindist, retentiontime, addch2toprec):
# method selection
	
	rtperm = False
	shufflefrag = False
	linmzperm = False
	sw_perm = False
	fragdb = False
		
	print("method used: ", method)
	
	if (method == "rtperm"):
		rtperm = True
		# saves already choosen rts
		rt_already = []
	elif (method == "shufflefrag"):
		shufflefrag = True
	elif (method == "linmzperm"):
		linmzperm = True
	elif (method == "sw_perm"):
		sw_perm = True
		# TODO: check if swath window files != empty
		sw = pd.read_csv(swath_windows, sep='\t')
	elif (method == "fragdb"):
		fragdb = True

# create temporary file	
	print('Create temporary file')
	# create validated temporary file 
	tmpfile = tempfile.NamedTemporaryFile(suffix='.tmp')
	fillTmpTSVWithValidTargetedExp(assay_library, tmpfile)
		
	# check output file extension
	base, extension = os.path.splitext(output)
	if extension == ".tsv" or extension == ".pqp":
		output_type = extension
	else:
		print('Output type is invalid! Please specify a valid output_type (.tsv, .pqp).')
		print('If output_type is not defined .pqp will be used automatically')
		output_type = ".pqp"
	
	out = base + output_type

# methods 
	# read tsv to dataframe
	df_assay_lib = pd.read_csv(tmpfile, sep = '\t')
	
	# add key "SW" to df
	if (sw_perm):
		df_assay_lib['SW'] = 0
		df_assay_lib['decoy_SW'] = 0
		df_assay_lib['new_mass_min'] = 0
		df_assay_lib['new_mass_max'] = 0

	unique_tgid = df_assay_lib.TransitionGroupId.unique()
	# divide df by unique TransitionGroupId
	d_lib = dict()
	d_lib_decoy = dict()
	for element in unique_tgid:
		d_lib[element] = df_assay_lib[df_assay_lib['TransitionGroupId'] == element]
		d_lib_decoy[element] = df_assay_lib[df_assay_lib['TransitionGroupId'] == element]

	if (fragdb): 
		fragdb = []
		fragdb = df_assay_lib['ProductMz'].tolist()
		fragdb.sort(reverse = True)

	# annotate swath window based on precursor mass and swath_windows file
	if (sw_perm):
		header = 1 # the file should be an header "start 	end"
		i = 0
		lowest_window = 0
		highest_window = sw.shape[0] - header
		for index, row in sw.iterrows():
			for key, value in d_lib_decoy.items():
				prec_mz = value['PrecursorMz'].unique()
				if (row[0] <= prec_mz[0] <= row[1]): 
					value['SW'] = i 
					current_window = i			
					# switch window to one lower or one higher - depending 
					# on precursor m/z value (nearer to the lower / upper window)
					if  (row[0] <= prec_mz[0] <= row[1] and current_window > lowest_window):
						new_window = current_window - 1
						value['decoy_SW'] = new_window
						mass_min = sw.iloc[new_window][0]
						mass_max = sw.iloc[new_window][1]
						value['new_mass_min'] = mass_min
						value['new_mass_max'] = mass_max - 1 # reduce overlap
					elif (current_window == lowest_window):
						new_window = current_window + 1
						#new_window = highest_window
						value['decoy_SW'] = new_window 
						mass_min = sw.iloc[new_window][0]
						mass_max = sw.iloc[new_window][1]
						value['new_mass_min'] = mass_min
						value['new_mass_max'] = mass_max -1 # no overlap
					else:
						print("There seems to be an error - please check the swath window file")
			i = i + 1

	# problem window size - if you have really big and really small windows - how to fit all
	# the deocys into the small windows? Approximately same amount of precursors - so it should not be
	# too bad.
	
	# add _decoy to TransitionGroupId, CompoundName
	for key, value in d_lib_decoy.items():
		value['TransitionGroupId'] = value['TransitionGroupId'].apply(str) + '_decoy'
		value['TransitionId'] = value['TransitionId'].apply(str) + '_decoy'
		value['CompoundName'] = value['CompoundName'].apply(str) + '_decoy'
		value['Decoy'] =  1
		value['Annotation'] = 'NA'

		if (fragdb):
			current_prec_mz = value['PrecursorMz'].unique()
			old_fragments = value['ProductMz'].tolist()
			new_fragments = old_fragments
			#print('current_precursor: ', current_prec_mz),
			#print('old: ',value['ProductMz'].tolist())
			# use all fragments lower than current precursor 
			number_frag = len(value['ProductMz'].tolist())
			# find index of closest element in list 
			closest_element = min(range(len(fragdb)), key=lambda i: abs(fragdb[i]-current_prec_mz))
			#print(closest_element)
			current_db = fragdb[closest_element:len(fragdb)]
			#print(current_db)
			while (new_fragments == old_fragments):
				new_fragments = random.sample(current_db, k=number_frag)
			value['ProductMz'] = new_fragments
			#print('new: ',value['ProductMz'].tolist())
					
		if (sw_perm):
			# randomly choose mz within a certain range  
			current_prec_mz = value['PrecursorMz'].unique()
			new_prec_mz = float(current_prec_mz)
			current_min = value['new_mass_min'].unique()
			current_max = value['new_mass_max'].unique()
			new_prec_mz = random.uniform(current_min, current_max)
			value['PrecursorMz'] = value['PrecursorMz'].replace([current_prec_mz], new_prec_mz)
			print("min: ", current_min, " ", "max: ", current_max)
			print("before :", current_prec_mz," ", "after: ", new_prec_mz)    

		if (rtperm):
			# retention time
			rt = value['NormalizedRetentionTime'].unique()
			newrt = float(rt)
			print("rt: ",rt)
			while (abs(newrt - rt) < rt_mindist): # loop as long as rt diff smaller than rt_mindist
				randomnum = random.random()
				print("rand:", randomnum)
				newrt = newrt + (randomnum - 0.5) * 400
				if (newrt in rt_already): # if it was used as decoy rt already
					newrt = rt
				if newrt < 0 or newrt > retentiontime: #has to be in chrom retention time range 
					newrt = rt # then it is again in the while loop
			value['NormalizedRetentionTime'] = value['NormalizedRetentionTime'].replace([rt], newrt)
			print("newrt: ",newrt)
			rt_already.append(newrt)

		if (shufflefrag):
			frag_mz = value['ProductMz'].tolist()
			print('frag_mz: ', frag_mz)
			newfrag_mz = frag_mz
			random.shuffle(newfrag_mz)
			print('newfrag_mz: ', newfrag_mz)
			value['ProductMz'] = newfrag_mz

		if (linmzperm):
			massdiffs = list()
			baseline = value['PrecursorMz'].unique()
			frag_mz = value['ProductMz'].tolist()
			frag_mz.sort(reverse = True)
			print('baseline:', baseline)
			print('frag_mz: ', frag_mz)
			for element in frag_mz:
				massdiff = round(float(baseline - element),6)
				massdiffs.append(massdiff)
				baseline = element
				
			# What happens if the same mass difference occures 
			# e.g. [27.01, 27.01] 
		
			newfraglist = frag_mz
			shuffled_massdiffs = list(massdiffs)
			
			# issue with "99.082647, 327.077562, 327.113947"
			# this step has to be performed on integer instead of float level
			# if all have the same mass difference add here an additional CH2 (to the first element)
			integer_massdiff = [ int(x) for x in shuffled_massdiffs ]
			unique = len(set(integer_massdiff)) == len(integer_massdiff) 
			if unique == False:
				print("INFO: Same mass difference occured")
				# find non unique value
				interger_massdiff_array = np.array(integer_massdiff)
				duplicate = [item for item, count in Counter(interger_massdiff_array).items() if count > 1]
				duplicate_index = integer_massdiff.index(duplicate[0])	
				# use one of the indices of duplicated items and add a CH2 
				shuffled_massdiffs[duplicate_index] = addCH2(shuffled_massdiffs[duplicate_index])
			
			# set random seed
			# random.seed(1)

			while(len(set(newfraglist).intersection(frag_mz)) > 1):
				random.shuffle(shuffled_massdiffs)
				print('shuffled_massdiff: ', shuffled_massdiffs)
				newfraglist = list() # reset in while loop
				base = value['PrecursorMz'].unique()
				for dmass in shuffled_massdiffs:
					newfrag = round(float(base - dmass),6)
					newfraglist.append(newfrag)
					base = newfrag
				if len(set(newfraglist).intersection(frag_mz)) == 1:
					element = list(set(newfraglist).intersection(frag_mz))[0]
					# get index of element in newfraglist 
					element_index = newfraglist.index(element)
					# add CH2 to element which is the same 
					newfraglist[element_index] = addCH2(newfraglist[element_index])
						
			# add CH2 additional to las element
			# newfraglist[-1] = addCH2(newfraglist[-1])

			print("target: ", frag_mz)
			print("decoy: ", newfraglist)
			
			# check if same number of entries
			if (len(frag_mz) != len(newfraglist)):
				print("Not same lenght")
				break

			value['ProductMz'] = newfraglist

		d_lib_decoy[key] = value
		
	# This works only for methods which only use the MS2 level
	# linmzperm, fragDB, shufflefrag
	
	print("addch2toprec: ", addch2toprec)
	
	if addch2toprec and method != "rtperm" and method != "sw_perm":
		for key, value in d_lib_decoy.items():
			value['PrecursorMz'] = addCH2(value['PrecursorMz'])
	

	# create df from dicts
	print("create dataframe") 
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

	print("Done")

if __name__ == "__main__":
    main()
    
