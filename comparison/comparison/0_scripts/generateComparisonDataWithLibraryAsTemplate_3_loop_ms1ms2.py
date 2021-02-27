# author: Oliver Alka
# date: 07.08.2019
#
# Script to generate a data for the quantification comparison of skyline and OpenSWATH/pyoprophet output 
# based on the compounds in the assay library used. 
#
# MS1_osw: aggr_prec_Peak_Area
# MS2_osw: Intensity
# MS1_sky: Area (precursor)
# MS2_sky: Area (transitions)


# packages
import os
import re
import math
import click
import statistics
import pandas as pd
import csv

allowed_rt_diff = 5

def rtMinToSec(x):
    return x*60

def extractBasename(x):
    base = os.path.basename(x)
    base = os.path.splitext(base)[0]
    return base

def doubleToInt(x):
    return int(x)
    
dilution_pattern = "Step(\\d+)"
def extractDilution(pattern,basename):
	return re.findall(pattern,basename)[0]
	
replicate_pattern = "Plasma(\\d)"
def extractReplicate(pattern,basename):
	return re.findall(pattern,basename)[0]
	
class MetaboliteResult:
	basename_sky = ""
	basename_osw = ""
	compounname = ""
	adduct = ""
	step = ""
	replicate = ""
	integer_mz = 0
	intensity_sky = 0.0
	intensity_ms1_sky = 0.0
	retentiontime_sky = float('NaN')
	intensity_osw = 0.0
	intensity_ms1_osw = 0.0
	retentiontime_osw = float('NaN')
	decoy = 0
	confusion = ""
	rt_diff = ""

number_of_dilutions = 10
number_of_replicates = 3

# command line tool options
@click.command()
@click.option('--input_template_assaylib', '-in_assay', envvar = 'input_template_assaylib', multiple = False, type = click.Path(), help = 'Assay library template which comparison is based on')
@click.option('--input_skyline', '-in_sky', envvar = 'input_skyline', multiple = False, type = click.Path(), help = 'Input of the skyline analysis (TransitionsResults)')
@click.option('--filtered_dir', '-filtered_dir', envvar = 'filtered_dir', multiple = False, type = click.Path(), help = 'Directory with qvalues filtered pyprophet results')
@click.option('--out_dir', '-out_dir', envvar = 'out_dir', multiple = False, type = click.Path(), help = 'Output directory for the comparison')


def main(input_template_assaylib,input_skyline,filtered_dir,out_dir):

	input_template_assaylib =input_template_assaylib
	input_skyline = input_skyline
	filtered_dir = filtered_dir
	outdir = out_dir

	for filename in os.listdir(filtered_dir):
		if filename.endswith(".tsv"):

			input_pyprophet = filtered_dir + filename
			output_tsv = outdir + "comp_" + filename
			print("input: ", input_pyprophet)
			print("output: ", output_tsv)

			assaylib = pd.read_csv(input_template_assaylib, sep = "\t")
			skyline = pd.read_csv(input_skyline, sep = ",")
			pyprophet = pd.read_csv(input_pyprophet, sep = "\t")

			# dataframe preparation
			skyline['rt_sec'] = skyline['RetentionTime'].map(rtMinToSec) 
			pyprophet['basename'] = pyprophet['filename'].map(extractBasename)
			skyline['integer_mz'] = skyline['PrecursorMz'].map(doubleToInt)
			pyprophet['integer_mz'] = pyprophet['mz'].map(doubleToInt)
			assaylib['integer_mz'] = assaylib['PrecursorMz'].map(doubleToInt)

			# key = compoundname_adduct_integermz_stepX_repX
			# dictionary(key, value(MetaboliteResult))
			entry_exists = []
			metabolites = {}
			for index, row in assaylib.iterrows():
				current_entry = row['CompoundName'] + '_' + str(row['integer_mz'])
				if current_entry not in entry_exists:
					entry_exists.append(current_entry)
					# generate_key and fill data structure
					for i in range(1,number_of_dilutions+1):
						for j in range(1,number_of_replicates+1):
							key = current_entry + "_" + "step" + str(i) + "_rep" + str(j)
							value = MetaboliteResult()
							value.compoundname = row['CompoundName']
							value.adduct = row['Adducts']
							value.step = i
							value.replicate = j
							value.integer_mz = row['integer_mz']
							if row['CompoundName'].endswith("_decoy"):
								value.decoy = 1
							metabolites[key]=value

			# fill with skyline results 
			missing_keys_sky = []
			for index, row in skyline.iterrows():
				# parse basename in this case to get the step and rep
				current_step = extractDilution(dilution_pattern, row['Replicate'])
				current_rep = extractReplicate(replicate_pattern, row['Replicate'])
				skykey = row['Peptide'] + '_' + str(row['integer_mz']) + "_" + "step" + str(current_step) + "_rep" + str(current_rep)
				# sum up intensity or mean rt based on the different transitions for the same metabolites (individual value in skyline)). 
				if skykey in metabolites:
					if row['Area'] !=0 and not math.isnan(row['Area']):
						current_value_sky = metabolites[skykey]
						current_value_sky.basename_sky = row['Replicate']
						# check if MS1 or MS2 level ("precursor" in column fragmentation)
						if "precursor" in row['FragmentIon']:
							current_value_sky.intensity_ms1_sky = current_value_sky.intensity_ms1_sky + row['Area']
							#print("---")
							#print(current_value_sky.intensity_ms1_sky)
						else:
							current_value_sky.intensity_sky = current_value_sky.intensity_sky + row['Area']	
							#print("---")
							#print(current_value_sky.intensity_sky)
						#if current_value_sky.retentiontime_sky != 0.0:
						if row['rt_sec'] != 0.0:
							#print("-------------------------")
							#print("current_value: ", current_value_sky.retentiontime_sky)
							#print("row_value: ", row['rt_sec'])
							#print("nan_check: ", math.isnan(current_value_sky.retentiontime_sky))
							if not math.isnan(current_value_sky.retentiontime_sky):
								current_value_sky.retentiontime_sky = statistics.mean([current_value_sky.retentiontime_sky, row['rt_sec']])
								#print("calc_stat: ", current_value_sky.retentiontime_sky)
							else:
								current_value_sky.retentiontime_sky = statistics.mean([row['rt_sec']])
								#print("calc_stat: ", current_value_sky.retentiontime_sky)
				
						metabolites[skykey] = current_value_sky
				else:
					missing_keys_sky.append(skykey)
	
			# fill with OSW results
			missing_keys_osw = []
			for index, row in pyprophet.iterrows():
				current_step = extractDilution(dilution_pattern, row['basename'])
				current_rep = extractReplicate(replicate_pattern, row['basename'])
				pypkey = row['compound_name'] + '_' + str(row['integer_mz']) + "_" + "step" + str(current_step) + "_rep" + str(current_rep)
				if pypkey in metabolites:
					current_value_osw = metabolites[pypkey]
					current_value_osw.basename_osw = row['basename']
					current_value_osw.intensity_osw = row['Intensity']
					current_value_osw.intensity_ms1_osw = row['aggr_prec_Peak_Area']
					current_value_osw.retentiontime_osw = row['RT']
					metabolites[pypkey] = current_value_osw
				else:
					missing_keys_osw.append(pypkey)
	
			#built entry
			# TODO: how to build the "new" confusion matrix using ms1 and ms2 information (?)
			header = ["basename_sky","basename_osw","compoundname","adduct","step","rep","integer_mz","int_sky","int_ms1_sky","int_pyp","int_ms1_pyp","rt_sky","rt_pyp","decoy","confusion","rt_diff"]
			lines = [header]
			for key, value in metabolites.items():

				# build confusion matrix based on basenames
				if value.basename_sky == value.basename_osw:
					value.confusion = "TP"
				if value.basename_sky != "" and value.basename_osw == "":
					value.confusion = "FN"
				if value.basename_sky == "" and value.basename_osw != "":
					value.confusion = "FP"
				if value.basename_sky == "" and value.basename_osw == "":
					value.confusion = "TN"
	
				# compare based on retention time 
				if value.retentiontime_sky > 0.0 and value.retentiontime_osw > 0.0:
					abs_rt_diff = abs(value.retentiontime_sky - value.retentiontime_osw)
					if abs_rt_diff > allowed_rt_diff:
						value.rt_diff = "diff"
						# if rt is different - set value of the confusion matrix to FP
						value.confusion = "FP"
						
				# Only if you want to find all combinations e.g. nan && 0.0 
				#elif math.isnan(value.retentiontime_sky) and value.retentiontime_osw > 0.0:
				#	value.rt_diff = "diff"
				#elif value.retentiontime_sky > 0.0 and math.isnan(value.retentiontime_osw):
				#	value.rt_diff = "diff"

				line = [value.basename_sky,
						value.basename_osw,
						value.compoundname,
						value.adduct,
						value.step,
						value.replicate,
						value.integer_mz,
						value.intensity_sky,
						value.intensity_ms1_sky,
						value.intensity_osw,
						value.intensity_ms1_osw,
						value.retentiontime_sky,
						value.retentiontime_osw,
						value.decoy,
						value.confusion,
						value.rt_diff]
				lines.append(line)

			with open(output_tsv, 'w+') as f:
				writer = csv.writer(f, delimiter='\t')
				for line in lines:
					writer.writerow(line)
			f.close()
			#print("Keys missing in skyline: " + ', '.join(map(str, missing_keys_sky)))
			#print("Keys missing in osw : " + ', '.join(map(str, missing_keys_osw)))
			print("Done")

	return 0
	
if __name__ == "__main__":
	main()