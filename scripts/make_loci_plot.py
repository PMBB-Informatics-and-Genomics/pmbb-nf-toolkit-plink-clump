import pandas as pd
import argparse as ap
import manhattan_plot
from manhattan_plot import *
import os

def make_parser():
    parser = ap.ArgumentParser()

    parser.add_argument('--analysis', required=True, help='clump analysis nickname')
    parser.add_argument('--loci', required=True, help='csv file with locus coordinates')
    parser.add_argument('--sumstats', nargs='+', required=True, help='summary stats files used for the clumping procedure')
    parser.add_argument('--annot', required=False, help='annotation file if included')

    return parser

# Parse the arguments
args = make_parser().parse_args()
analysis = args.analysis
sumstats_files = args.sumstats
annot_file = args.annot
loci_file = args.loci
with_annot = annot_file is not None

# Read in and concatenate all summary stats files
dfs = []
for f in sumstats_files:
    dfs.append(pd.read_table(f))
df = pd.concat(dfs)
print(df)

# Write it so that we can load it with the plotting package later
df.to_csv('temp_sumstats.tsv', sep='\t', index=False)

# Read in the biofilter annotations
if with_annot:
    annot_df = pd.read_csv(annot_file)
    annot_df = annot_df.rename(columns={'Gene': 'ID'})
    annot_df = annot_df.set_index('Var_ID')
    print(annot_df)

# Read in the loci coordinates and rename columns
loci_df = pd.read_csv(loci_file).rename(columns={
    'MIN_POS': 'START', 'MAX_POS': 'END',
    'CHR': '#CHROM', 'BP': 'POS'
})

loci_df['ID'] = annot_df.loc[loci_df['Lead_SNP'], 'ID'].values
print(loci_df.columns)
print(loci_df)

unique_loci = []
for locus_id, subDF in loci_df.groupby('ID'):
    if len(subDF) > 1:
        subDF['ID'] = subDF['ID'] + '.' + [str(i) for i in range(1, len(subDF)+1)]
    unique_loci.append(subDF)

loci_df = pd.concat(unique_loci)
print(loci_df)

# Start the Manhattan plot object
mp = ManhattanPlot('temp_sumstats.tsv', title=f'Significant GWAS Loci for {analysis.replace("_", " ")}')
mp.load_data(delim='\t')
mp.clean_data(col_map={'CHR': '#CHROM', 'BP': 'POS', 'SNP': 'ID'})

# Add annotations if needed and thin
if with_annot:
    mp.add_annotations(annot_df, extra_cols=['RSID'])
mp.get_thinned_data()

# Set plotting parameters
mp.update_plotting_parameters(
    vertical=True,
    merge_genes=True
)

# Make a signals-only ploy
mp.signal_plot_with_specific(
    loci_df,
    extra_cols={'RSID': 'RSID'},
    save=f'{analysis}.loci.png'
)

os.remove('temp_sumstats.tsv')