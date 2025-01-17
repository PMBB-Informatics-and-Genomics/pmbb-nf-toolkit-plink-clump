import sys
import pandas as pd
import argparse as ap
import numpy as np
import os

def make_parser():
    parser = ap.ArgumentParser()

    parser.add_argument('--clumps', nargs='+', help='Output from the plink clump call')
    parser.add_argument('--extract-files', nargs='+', help='variant IDs used in the clumping procedure')
    parser.add_argument('--analysis', help='nickname of the clumping analysis')

    return parser

# Parse the arguments
args = make_parser().parse_args()
clump_files = args.clumps
analysis = args.analysis
extract_files = args.extract_files

# Build the output file names
output_clumps = f'{analysis}.clumps.csv'
output_loci = f'{analysis}.loci.csv'

# Iterate over clump files
dfs = []
print(clump_files)
for f in clump_files:
    if os.path.getsize(f) == 0:
        print(f)
        continue
    dfs.append(pd.read_table(f, sep='\s+'))

# Concatenate clump files together
if len(dfs) == 0:
    all_clumps = pd.DataFrame(columns=['#CHROM', 'POS', 'ID', 'MIN_POS', 'MAX_POS'])
else:
    all_clumps = pd.concat(dfs)
all_clumps = all_clumps.set_index('ID')
print(all_clumps)
print(all_clumps.columns)

# Read in list of variants used in clumping
extract_labels = pd.concat([pd.read_table(f, header=None, index_col=3) for f in extract_files])
print(extract_labels)

# This is where we add the minimum and maximum position
for lead_snp, row in all_clumps.iterrows():
    print(row)
    clump_vars = row['SP2'].replace('(1)', '').replace('NONE', '').replace('.', '').split(',')
    clump_vars.append(lead_snp)
    clump_vars = [var for var in clump_vars if var != '']
    print(clump_vars)

    clump_pos = extract_labels.loc[clump_vars]
    print(clump_pos)

    all_clumps.loc[lead_snp, 'MIN_POS'] = clump_pos[1].min().astype(int)
    all_clumps.loc[lead_snp, 'MAX_POS'] = clump_pos[1].max().astype(int)

# Write output of all clumps
print(all_clumps)
all_clumps = all_clumps.sort_values(by=['#CHROM', 'MIN_POS'])
all_clumps['ANALYSIS'] = analysis
all_clumps.index.name = 'Lead_SNP'
all_clumps.to_csv(output_clumps)

# Now merge physically overlapping clumps
all_loci = all_clumps.copy()
overlap_test = np.logical_and(all_loci['MIN_POS'].iloc[1:].values < all_loci['MAX_POS'].iloc[:-1].values, all_loci['#CHROM'].iloc[1:].values == all_loci['#CHROM'].iloc[:-1].values)

while np.any(overlap_test):
    subDF1 = all_loci.iloc[:-1][overlap_test]
    subDF2 = all_loci.iloc[1:][overlap_test]

    merge_row1 = subDF1.iloc[0]
    merge_row2 = subDF2.iloc[0]

    # Keep lead SNP information with the smaller P value
    winner_row = merge_row1 if merge_row1['P'] < merge_row2['P'] else merge_row2

    new_row = winner_row.copy()
    new_row['SP2'] = merge_row1['SP2'] + merge_row2['SP2']
    new_row['TOTAL'] = merge_row1['TOTAL'] + merge_row2['TOTAL']
    new_row['MIN_POS'] = min(merge_row1['MIN_POS'], merge_row2['MIN_POS'])
    new_row['MAX_POS'] = max(merge_row1['MAX_POS'], merge_row2['MAX_POS'])

    all_loci = all_loci.drop([merge_row1.name, merge_row2.name])
    all_loci.loc[winner_row.name] = new_row
    all_loci = all_loci.sort_values(by=['#CHROM', 'MIN_POS'])
    overlap_test = np.logical_and(
        all_loci['MIN_POS'].iloc[1:].values < all_loci['MAX_POS'].iloc[:-1].values,
        all_loci['#CHROM'].iloc[1:].values == all_loci['#CHROM'].iloc[:-1].values)

print(all_loci)
all_loci.to_csv(output_loci)
