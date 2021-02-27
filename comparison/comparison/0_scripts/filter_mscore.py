#
# author: Oliver Alka
# date: 20190013
#
# Filter pyprophet output at differenct m / q-values
#

import os
import click 
import pandas as pd
import numpy as np

def extractBasename(x):
    base = os.path.basename(x)
    base = os.path.splitext(base)[0]
    return base

# command line tool options
@click.command()
@click.option('--input_pyprophet', '-in', envvar = 'input_pyprophet', multiple = False, type = click.Path(), help = 'Unfiltered pyprophet results')
@click.option('--outdir', '-outdir', envvar = 'outdir', multiple = False, type = click.Path(), help = 'Output directory for filtered results')

def main(input_pyprophet,outdir):

	input_unfiltered_pyprophet = input_pyprophet
	out_dir = outdir

	# 30% FDR to 1% FDR
	qvalues = [float(np_float) for np_float in np.flip((np.arange(0.01,0.31,0.01)))]
	
	# append for smaller gaps in the filtering
	qvalues = qvalues + [float(np_float) for np_float in np.flip((np.arange(0.0001,0.01,0.001)))]
	print(qvalues)
	
	basename = extractBasename(input_unfiltered_pyprophet)

	unfiltered = pd.read_csv(input_unfiltered_pyprophet, sep = "\t")
	
	# pre-filter by peak_group_rank == 1 
	prefiltered = unfiltered[unfiltered['peak_group_rank'] == 1]

	for element in qvalues:
		filtered = prefiltered
		filtered = filtered[filtered['m_score'] < element]

		name = basename + "_" + str(element).replace(".","")+ ".tsv"
		filtered.to_csv(out_dir+name,sep = "\t")
		
	

	print("Done")
	
	return 0
	
if __name__ == "__main__":
	main()

