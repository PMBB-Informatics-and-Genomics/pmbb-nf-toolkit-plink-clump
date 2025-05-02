import sys
import pandas as pd
import argparse as ap
import numpy as np
import os

def make_parser():
    parser = ap.ArgumentParser()
    parser.add_argument('--clumps', help='Output from the plink clump call')

    return parser

# Parse the arguments
args = make_parser().parse_args()
clump_file = args.clumps

# Read the clump file and write the lead snps
output_snps = 'clump_lead_snps.txt'
df = pd.read_table(clump_file, sep='\s+')
open(output_snps, 'w+').write('\n'.join(df['ID'].unique()) + '\n')

