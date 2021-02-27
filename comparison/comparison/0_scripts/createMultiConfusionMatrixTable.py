#
# author: Oliver Alka
# date: 13.09.2019 
#
# Script to create table from multiple thresholds / FDRs 
# pyprophet results wihtout FDR have to be filtered using the "filter_mscore.py" 
# script

import os
import re
import click
import pandas as pd 

#workaround to insert "." in estimated FDR after importing 
def insert_str(string, str_to_insert, index):
    return string[:index] + str_to_insert + string[index:]

def extractBasename(x):
    base = os.path.basename(x)
    base = os.path.splitext(base)[0]
    return base

eFDR_pattern = "\d+"
def extractEstimatedFDR(pattern,basename):
    eFDR = "NA"
    eFDR = re.findall(pattern,basename)
    return eFDR
    
# command line tool options
@click.command()
@click.option('--input_dir', '-input_dir', envvar = 'input_dir', multiple = False, type = click.Path(), help = 'Directory with compared skyline and pyprophet results')
@click.option('--out', '-out', envvar = 'out', multiple = False, type = click.Path(), help = 'Output of confusion matrix (.tsv)')

def main(input_dir,out):

	input_dir = input_dir
	output = out

	# create df to store the stuff 
	confusionmatrixes = pd.DataFrame(columns=['Name', 'TP', 'FP', 'TN', 'FN', 'estimatedFDR'])
	print(confusionmatrixes)

	for filename in os.listdir(input_dir):
		if filename.endswith(".tsv"):
			basename = extractBasename(filename)
			estimatedFDR = str(extractEstimatedFDR(eFDR_pattern,basename)[-1])
			estimatedFDR = insert_str(estimatedFDR,".",1)
			print(basename)
			print(estimatedFDR)

			# read df
			current_file = pd.read_csv(input_dir+filename, sep = "\t")
			# remove decoys
			current_file = current_file[current_file['decoy'] == 0]
			# count occurance and add to confusionmatrixes
			countocc = current_file['confusion'].value_counts(dropna=False)        
			print(countocc) 
			print(type(countocc))
			
			# correct missing values 
			if 'TP' in countocc:
				tp = countocc['TP']
			else:
				tp = 0
				
			if 'FP' in countocc:
				fp = countocc['FP']
			else:
				fp = 0
				
			if 'TN' in countocc:
				tn = countocc['TN']
			else:
				tn = 0
				
			if 'FN' in countocc:
				fn = countocc['FN']
			else:
				fn = 0
			
			dfrow = ({'Name':basename,
					  'TP':tp,
					  'FP':fp,
					  'TN':tn,
					  'FN':fn,
					  'estimatedFDR':estimatedFDR})
			#print(dfrow)
			confusionmatrixes = confusionmatrixes.append(dfrow, ignore_index=True)

	confusionmatrixes.to_csv(output,sep="\t",index=False)


	print("Done")
	
	return 0

if __name__ == "__main__":
	main()