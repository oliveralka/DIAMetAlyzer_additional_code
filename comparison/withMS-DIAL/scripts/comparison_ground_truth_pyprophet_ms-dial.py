# author: Oliver Alka
# date: 28.06.2021
#
# Script to compare the fitered MS-DIAl peak list results to the ground truth and pyprohet (e.g. 5% FDR) results.
#

# packages
import os
import click
import pandas as pd

def format_adduct(adduct):
    m_adduct = adduct[0:-1]
    m_charge = adduct[-1]
    adduct = "[" + m_adduct + "]" + m_charge
    return adduct

# command line tool options
@click.command()
@click.option('--input_gd', '-input_gd', envvar = 'input_gd', multiple = False, type = click.Path(), help = 'Input of ground truth and pyprophet confusion matrix.')
@click.option('--input_md', '-input_md', envvar = 'input_md', multiple = False, type = click.Path(), help = 'Input of filtered MS-DIAL results')
@click.option('--output_comparison', '-output_comparison', envvar = 'output_comparison', multiple = False, type = click.Path(), help = 'Output of comparison between the ground truth, pyprophet and MS-DIAL')
@click.option('--allowed_rt_diff', '-allowed_rt_diff', envvar = 'allowed_rt_diff', multiple = False, type = float, default = 5.0, help = 'Allowed retention time deviation. If the retention time defiation is too large the compound will be regarded as false positive (FP).')

def main(input_gd,input_md,output_comparison,allowed_rt_diff):

    # extract basename based on ms-dial file input
    basename_md = "_".join(os.path.splitext(input_md)[0].split("_")[-2:]).split("/")[1]
    print(basename_md)

    # read inputs
    gd = pd.read_csv(input_gd, sep = "\t")
    md = pd.read_csv(input_md, sep = "\t")

    # add basename to md 
    md['basename_md'] = basename_md

    # filter rows with empty basename or wrong basename (e.g. from one of the other dilutions)
    gd.dropna(0,subset=['basename_sky','basename_osw'],inplace=True)
    for index, row in gd.iterrows():
        if row['basename_sky'] != basename_md:
            gd.drop(index, inplace = True)

    # reformat adduct to [M+x]+ format for easier comparison 
    gd['adduct'] = gd['adduct'].map(format_adduct)

    md['key'] = md['Title'] + "_" + md['Adduct']
    gd['key'] = gd['compoundname'] + "_" + gd['adduct']

    # merge both dataframes into new "mer" daraframe
    mer = pd.merge(gd, md, how='left', on='key')

    # initialize additional columns
    mer['confusion_md'] = ""
    mer['rt_diff_md'] = ""

    # TP, FN, FP, TN based on basename (similar to before)
    for index, row in mer.iterrows():
        if row['basename_sky'] == row['basename_md']:
            mer.at[index, 'confusion_md'] = "TP"
        else:
            mer.at[index, 'confusion_md'] = "FN"

        abs_rt_diff = abs(row['rt_sky'] - row['RT (sec)'])
        if abs_rt_diff > allowed_rt_diff:
            mer.at[index, 'rt_diff_md'] = "diff"
            mer.at[index, 'confusion_md'] = "FP"

    mer.to_csv(output_comparison, sep = "\t")

    return 0

if __name__ == "__main__":
    main()
