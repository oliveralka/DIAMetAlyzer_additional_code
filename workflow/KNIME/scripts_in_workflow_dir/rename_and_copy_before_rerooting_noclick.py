# author: Oliver Alka
# date: 25.05.2019
#
# read *.dot files in dir, rename with ScanNr and Name form path
# and copy them to different location
# for passatuto processing

import os 
import fnmatch 
import re
import shutil

def find(pattern, path):
	result = []
	for root, dirs, files in os.walk(path):
		for name in files:
			if fnmatch.fnmatch(name, pattern):
				result.append(os.path.join(root, name))
	return result

def main(input,output):
	sep = os.path.sep
	treepaths = find('*.dot', input)
	names = []
	for path in treepaths:
		# split path first
		splitpath = os.path.split(path)[0]
		# split by path separator name should always be in second last element
		splitpath = splitpath.split(sep)[-2]
		# extract compound_name and scanNr - append it after 1Here_C4H10NO3PS_M+H+
		# pattern -/
		splitpath = splitpath.rsplit('-', 2)
		posname = splitpath[-2]
		name = splitpath[-1]

		# posname can either be a string of numbers, or parts of the name
		# 'CycloateRO', 'NEET1915'
		# '676005', 'Fenamiphos1680' 

		try:
			int(posname) 
		except ValueError:
			name = posname + "-" + name

		names.append(name)

	if len(names) == len(treepaths):
		for i in range(len(treepaths)):
			basename = os.path.basename(treepaths[i])
			# add name on second positions (after 1))
			newname = basename[:1] + names[i] + basename[1:]
			# copy the renamed file to the specific directory
			shutil.copy(treepaths[i],output+newname)
	else:
		print("The two lists have not the same lenght, there seems to be something wrong")
		
	
		
