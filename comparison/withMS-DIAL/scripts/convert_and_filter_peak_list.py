# author: Oliver Alka
# date: 28.06.2021
#
# Script to filter the MS-DIAl peak list to compounds with MS2 references and
# convert the retention time to seconds
#

# packages
import os
import click
import pandas as pd

# functions
def rtMinToSec(x):
    return x*60

# command line tool options
@click.command()
@click.option('--input_dir', '-input_dir', envvar = 'input_dir', multiple = False, type = click.Path(), help = 'Directory with MS-DIAL peak list exported results')
@click.option('--output_dir', '-output_dir', envvar = 'output_dir', multiple = False, type = click.Path(), help = 'Output directory for the filtered and converted MS-DIAL peak list results')

def main(input_dir, output_dir):

    input_dir = input_dir
    output_dir = output_dir

    for filename in os.listdir(input_dir):
        if filename.endswith(".txt"):

            input_file = input_dir + filename
            print("Input: "  +input_file)

            peaklist = pd.read_csv(input_file, sep = "\t")

            # convert rt to seconds and replace the current min columns
            peaklist['RT (min)'] = peaklist['RT (min)'].map(rtMinToSec)
            peaklist['RT left(min)'] = peaklist['RT left(min)'].map(rtMinToSec)
            peaklist['RT right (min)'] = peaklist['RT right (min)'].map(rtMinToSec)
            peaklist['Reference RT'] = peaklist['Reference RT'].map(rtMinToSec)

            peaklist.rename(columns={"RT (min)":"RT (sec)", "RT left(min)":"RT left(sec)", "RT right (min)":"RT right (sec)"}, errors="raise", inplace=True)

            # filter based on reference (remove "Unknown" and "w/o ...")
            for index, row in peaklist.iterrows():
                if row['Title'] == 'Unknown' or row['Title'].startswith('w/o'):
                    peaklist = peaklist.drop(index)

            current_output_file = output_dir + filename
            if os.path.exists(current_output_file):
                print("Error: The current output file " + filename + " does already exist. Please specify a different output directory, or delete the file.")
            else:
                print("Output: " + current_output_file)
                peaklist.to_csv(current_output_file, sep = "\t")

    return 0

if __name__ == "__main__":
    main()




