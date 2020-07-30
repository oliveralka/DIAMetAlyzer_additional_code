#Pesticides_UIS

#Author: Premy Shanthamoorthy
#Contact: premy.shan@gmail.com
#Date: 07/29/2020

"""
Evaluating complex backgrounds that may cause ambiguities in the measurement
of metabolites. A list of pesticides are compared against the NIST Tandem MS/MS 
library. This script profiles different methods for unique transitions as follows
in the Q1 and Q3 phases with the mentioned filters to identify the number of unique
ion signatures (UIS) per mol_id.

Q1/Q3 used
MS1 -   25 ppm / -
MRM -   0.7 Da / 0.7 Da 
SWATH - Variable Windows / 25 ppm

This tool will also measure the number of interferences for each transition
(the number of identical transitions within the range of metabolites filtered
as specified above). 
"""

import pandas as pd
import numpy as np
import heapq
import rdkit
import difflib
import re
import itertools
import time
import matplotlib.pyplot as plt
import seaborn as sns
import pyper as pr
from math import sqrt
import scipy.spatial as sp
from statistics import mode
import ast

"""
function read:
Reads the NIST dataset compounds and spectra

input: .pkl files (a list of compounds, a list of spectra)
output: pandas dataframes of compounds and spectra (cf and spectra)
"""
def read(compounds, spectra):
    allcomp = pd.read_pickle(compounds) 
    spectra = pd.read_pickle(spectra) 
    cf = allcomp.dropna(subset = ['inchi', 'smiles', 'cas_num']) 
    compounds_null = allcomp.loc[allcomp['inchi'].isnull()]
    return cf, spectra

"""
function read_pes:
Reads in full pesticides dataset and removes pesticides that have <3 transitions

input: .csv file (full pesticides dataset)
output: pandas dataframes of pesticides
"""
def read_pes(filename):
    pesticides = pd.read_csv(filename)
    pesticides = pesticides[['CompoundName','PrecursorMz', 'ProductMz', 'LibraryIntensity']]
    pesticides['tuple'] = list(zip(pesticides['ProductMz'], pesticides['LibraryIntensity']))
    pesticides = pesticides.groupby(['PrecursorMz','CompoundName'])['tuple'].apply(list).reset_index(name='peaks')
    pesticides.columns = ['prec_mz', 'CompoundName','peaks']

    subtract = []
    
    for i, row in pesticides.iterrows():
        if (pesticides['peaks'][i].__len__()<3): #only keeping pesticides with >= 3 transitions
            subtract.append(i)
                    
    for x in subtract:
        pesticides = pesticides.drop(x)

    return(pesticides)

"""
function filter:
Filters the compound list based on the inchikey to take out chiral isomers but keep structural isomers. 

input: List of compounds from NIST
output: Filtered compound list
"""

def filter_comp(compounds_filt):
    subtract = [] #ones that have chiral isomers that need to be removed #6378
    compounds_filt = compounds_filt.drop_duplicates(subset='inchi')
    
    for index, row in compounds_filt.iterrows():
        sample_set = compounds_filt.loc[compounds_filt['exact_mass'] == row['exact_mass']] # specific isomers of the compound being looked at
        sample_set.drop(index)
        inchi = row['inchi']
        matches = []

        if len(sample_set)>1: #there are isomers
            if '/m' in inchi:
                connectivity = ((re.search(('InChI(.*)/m'), inchi)).group(1))
            else:
                connectivity = ""
                                            
            for i2, isomer in sample_set.iterrows(): #check if stereo connectivity is the same, if so, a chiral isomer so can be removed
                other = isomer['inchi']
                if '/m' in other:
                    check = ((re.search(('InChI(.*)/m'), other)).group(1))
                else:
                    check = other
                if (check == connectivity) and (i2 not in subtract):
                    subtract.append(i2)
                    matches.append(other)

    for x in subtract:
        compounds_filt = compounds_filt.drop(x)
    
    return compounds_filt

"""
function choose_background_and_query:
Per mol_id, choosing a specific query (based on given mol_id and parameters), followed by choosing the background
to compare this query to based on the Q1 and Q3 filters.

Input:

mol_id:           index of pesticide
molid_same:       molecular id of pesticides that are also measured in NIST 
Q1 parameters:    change, ppm - if ppm is filled, that will take priority over change
query parameters: col_energy, col_gas, ion_mode, inst_type, adduct
Q3 parameters:    If q3 = True, Q3 parameters will be taken into account - change_q3 or ppm_q3 parameters (if ppm_q3 is filled, that will take priority over change_q3)
top_n:            the top % of fragment ions to look at from which the most intense is chosen for the transition) 
variable:         Boolean, True if variable windows are used (if True, change and ppm are ignored)

Output: query (row), background_filt (background for query prec_mz), transitions_q1 (background for frag_mz),
        query_frag_mz_value (value of transition frag_mz), query_frag_mz (query fragment m/z with intensity)
"""

def choose_background_and_query_pesticides(spectra_filt, background2, mol_id, molid_same,change=0.7, ppm=0, change_q3=0.7, ppm_q3 = 0, q3= False, top_n= 0.1, UIS_num = 1, plot_back =False, variable=False):
    query = spectra_filt.loc[(spectra_filt.index == mol_id)]

    background_filt = spectra_filt.drop(query.index) #drop spectra from NIST if same pesticide is measured
    if molid_same != -1:
        same = background2.loc[background2['mol_id']==molid_same.values[0]]
        background2 = background2.drop(same.index)

    inst_type = ['IT/ion trap', 'Q-TOF', 'HCD','IT-FT/ion trap with FTMS']
    background2 = background2.loc[background2['inst_type'].isin(inst_type)]
    background2 = background2[['prec_mz', 'peaks']]
    background_filt = pd.concat([background_filt, background2])

    if len(query) == 1:
        #choosing a transition (highest intensity for that query)
        query_prec_mz = float(query['prec_mz'].item())

        #choosing background
        #for variable windows
        if variable==True:
            if (query_prec_mz >= 99.5) & (query_prec_mz <= 157.7):
                low = 99.5
                high = 157.7
            elif (query_prec_mz >= 156.7) & (query_prec_mz <= 242.6):
                low = 156.7
                high = 242.6
            elif (query_prec_mz >= 241.6) & (query_prec_mz <= 370.3):
                low = 241.6
                high = 370.3
            elif (query_prec_mz >= 369.3) & (query_prec_mz <= 465.9):
                low = 369.3
                high = 465.9
            elif (query_prec_mz >= 464.9) & (query_prec_mz <= 511.7):
                low = 464.9
                high = 511.7
            elif (query_prec_mz >= 510.7) & (query_prec_mz <= 557.5):
                low = 510.7
                high = 557.5
            elif (query_prec_mz >= 556.5) & (query_prec_mz <= 627.7):
                low = 556.5
                high = 627.7
            elif (query_prec_mz >= 626.7) & (query_prec_mz <= 886.1):
                low = 626.7
                high = 886.1            
        else:          
            if ppm != 0:
                change = ppm/1000000
            low = query_prec_mz - change/2
            high = query_prec_mz + change/2
        background_filt = background_filt.loc[background_filt['prec_mz'].between(low, high, inclusive = False)]

        #choosing the fragment
        query_frag_mz =  list(query['peaks'])[0]
        query_frag_mz.sort(key = lambda x: x[1], reverse = True)
    
        if query_frag_mz[0][0] != int(query_prec_mz):
            start = 0
        else:
            start = 1
            UIS_num += 1

        query_frag_mz = query_frag_mz[start:UIS_num]
        query_frag_mz_values = [query[0] for query in query_frag_mz]

        if q3 == True:
            for transition in query_frag_mz_values:
                if ppm_q3 != 0:
                    change_q3 = ppm_q3/1000000.0
                low = transition - change_q3/2.0
                high = transition + change_q3/2.0

                transitions_q1 = [[(a,b) for (a,b) in peaklist if a>low and a<high and (b>(1*top_n))] for peaklist in background_filt['peaks']] #do transitions here    
                transitions_q1 = [x for x in transitions_q1 if x!= []]
                transitions_q1 = list(itertools.chain.from_iterable(transitions_q1))
                transitions_q1.sort(key = lambda x: x[1], reverse = True)

                background_filt = background_filt.loc[(background_filt['peaks'].apply(lambda x: any(transition in x for transition in transitions_q1)))]
        
    else:
        query_frag_mz = 0
        query_frag_mz_values = 0

    if plot_back == True:
        ax = sns.distplot(list(background_filt['prec_mz']),bins = 20,kde = False)
        ax.grid()
        plt.show()
    
    return query, background_filt, query_frag_mz_values, query_frag_mz  
    

"""
function profile:
Based on the given parameters calculates the number of UIS and Interferences by mol_id.

Input: parameters for choose_background_and_query_pesticides
Output: pandas dataframe of compounds list with added columns of 'UIS' and 'Interferences'
"""

def profile_pesticides(compounds_filt, spectra_filt, background2, change = 0, ppm = 0, change_q3 = 0, ppm_q3 = 0, q3 = True, top_n = 0.1,
            mol_id = 0, plot_back = False, UIS_num=1, variable=False):
    UIS = []
    Interferences = []
    copy = spectra_filt
    
    for i, molecule in spectra_filt.iterrows():
        molid = i
        name = molecule['CompoundName']
        molid_same = compounds_filt.loc[compounds_filt['name']==name]
        if molid_same.empty == True:
            molid_same=-1
        else:
            molid_same=molid_same.index
        query, background, frag_mz, frag_mz_int = choose_background_and_query_pesticides(mol_id = molid, molid_same = molid_same,change = change, background2 = background2, ppm = ppm, change_q3 = change_q3, ppm_q3 = ppm_q3, q3 = q3, top_n = top_n, spectra_filt = spectra_filt, UIS_num=UIS_num, variable=False)
        if len(query) != 0:
            interferences = len(np.unique(background.index))
            if interferences == 0:
                unique = 1
            else:
                unique = 0
                
            UIS.append(unique)
            Interferences.append(interferences)
        else:
            UIS.append(-1) #if no query was found
            Interferences.append(-1)
    copy['UIS'] = UIS
    copy['Interferences'] = Interferences
    return copy


"""
function pesticides_profiler:
Profiles datasets according to Q1/Q3 window sizes:

Input:
change    - MS1 window size in Da, 0 if not in Da
change_q3 - MS2 window size in Da, 0 if not in Da 
ppm       - MS1 window size in ppm, 0 if not in ppm
ppm_q3    - MS2 window size in ppm, 0 if not in ppm
q3        - Boolean Value, True if Q3 window is used, False if not 
filename  - dataset that you would like to compare in .csv formt (i.e. pesticides)
compounds - list of compounds from additional dataset to compare, in .pkl format (i.e. NIST)
spectra   - list of spectra from additional dataset to compare, in .pkl format (i.e. NIST)

Output: Pandas dataframe with original dataset (filename) with added columns of 'UIS' and 'Interferences'
"""

def pesticides_profiler(filename, compounds, spectra, change=25, ppm=0, change_q3 =0, ppm_q3=0, top_n=0.1, mol_id=0, UIS_num=1, q3=False, variable=False):
    pesticides = read_pes(filename=filename)
    allcomp, spectra = read(compounds = compounds, spectra = spectra)#NIST data
    compounds_filt_NIST = filter_comp(compounds_filt = allcomp)
    spectra_filt_NIST = spectra.loc[spectra['mol_id'].isin(list(compounds_filt_NIST.mol_id))]
    spectra_filt_NIST.peaks = [[(a,b/1000) for (a,b) in peaklist] for peaklist in spectra_filt_NIST['peaks']]

    start = time.time()
    profiled = profile_pesticides(compounds_filt = compounds_filt_NIST, change = change, ppm = ppm, change_q3 = change_q3, ppm_q3 = ppm_q3, q3 = q3, background2 = spectra_filt_NIST, top_n = 0.1, mol_id = mol_id, spectra_filt = pesticides, UIS_num = UIS_num, variable=variable)
    profiled_filtered = profiled.loc[profiled['Interferences'] != -1]
    end = time.time()
    list_mol_ids = list(profiled_filtered.index)
    print("The unique identities and interferences for all mol_id will now be shown for profiling")
    print("The number of unique mol_id is: " + str(len([x for x in profiled['UIS'] if x == 1])))
    print("Time to completion of profiler: " + str(end-start))
    return profiled




