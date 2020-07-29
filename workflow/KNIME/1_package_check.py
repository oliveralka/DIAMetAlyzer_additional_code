import importlib.util
import sys 

package_names = ['numpy',
'pandas',
'pyprophet',
'pyopenms'
]

for package in package_names:
	spec = importlib.util.find_spec(package)
	if spec is None:
		print(package +" is not installed")
		print("Please insall the package using:") 
		print("pip install " + package)

# hack for the node to finish without an error on the KNIME side
import pandas
output_table = pandas.DataFrame({'col1': [1, 2], 'col2': [3, 4]})


############
# All packages used within the project:
# importlib, sys, os, pandas, tempfile, subprocess, 
# io, glob, re, numpy, pyopenms, csv, collections, 
# fnmatch, shutil
