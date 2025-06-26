import sys
import pandas as pd
import argparse as ap
import numpy as np

def make_parser():
    parser = ap.ArgumentParser()

    parser.add_argument('--sumstats', help='Summary stats file')
    parser.add_argument('--plink-bim', help='Chr-specific bim file')
    parser.add_argument('--pfile', action='store_true', help='plink flag for pgen file set')
    parser.add_argument('--bfile', action='store_true', help='plink flag for bed/bim/fam file set')
    parser.add_argument('--analysis', help='analysis nickname')
    parser.add_argument('--chr', help='chromosome')
    parser.add_argument('--p1thresh', type=float, help='clump p1 (primary threshold)')
    parser.add_argument('--p2thresh', type=float, help='clump p2 (secondary threshold)')
    parser.add_argument('--bp-window', type=int, help='bp window for searching around lead snps (ex 500000)')
    parser.add_argument('--kb-window', type=int, help='bp window for searching around lead snps (ex 500000)')

    parser.add_argument('--chr-col', help='chromosome column')
    parser.add_argument('--pos-col', help='position column')
    parser.add_argument('--id-col', help='variant ID column')
    parser.add_argument('--p-col', help='p-value column')
    parser.add_argument('--a1-col', help='A1 column')
    parser.add_argument('--a2-col', help='A2 column')

    return parser

# Parse arguments
args = make_parser().parse_args()
chr = args.chr
plink_bim = args.plink_bim

is_pfile = args.pfile and not args.bfile

analysis = args.analysis
sumstats = args.sumstats
p1thresh = args.p1thresh
p2thresh = args.p2thresh
bp_window = args.bp_window if args.bp_window is not None else 1000*args.kb_window

# Build output file names
output_sumstats = f'{analysis}.{chr}.txt.gz'
output_extract = f'{analysis}.{chr}.extract.txt'
output_min_p = f'{analysis}.{chr}.min_p.csv'

# Read in the bim file
if not is_pfile:
    bim = pd.read_table(plink_bim, sep='\s+', header=None, dtype=str)
else:
    bim = pd.read_table(plink_bim, sep='\s+', header=None, dtype=str, comment='#')
    pvar_to_bim_col_map = {
        0: 0, # CHROM -> CHR
        1: 3, # POS -> BP
        2: 1, # ID -> SNP
        3: 5, # REF -> A2
        4: 4, # ALT -> A1
        5: 'Other'
    }
    bim = bim.rename(columns=pvar_to_bim_col_map)

print(bim.head())

# Set up two versions for allele directionality
bim1 = bim.set_index([0, 3, 4, 5])
bim2 = bim.set_index([0, 3, 5, 4])

col_map = {
    args.chr_col : 'CHR',
    args.pos_col : 'BP',
    args.id_col : 'VAR_ID',
    args.a1_col : 'A1',
    args.a2_col : 'A2',
    args.p_col : 'P'
}

# Read in the summary stats file in chunks because it's big
ss_dfs = []
for df in pd.read_table(sumstats, sep='\s+', chunksize=1E6, dtype={args.chr_col: str, args.pos_col: str}):
    df = df.rename(columns=col_map)
    df['CHR'] = df['CHR'].astype(str).str.replace('chr', '')
    df = df[df['CHR'].astype(str) == str(chr)]
    if len(df) == 0:
        continue

    # Filter on p value
    df = df[df['P'] <= p2thresh]

    # Set the index
    df = df.set_index(['CHR', 'BP', 'A1', 'A2'])
    print(df.head())

    # Match on chromosome, position, and alleles
    forward_match = df.index.intersection(bim1.index)
    backward_match = df.index.intersection(bim2.index)

    print('Forward Match:', len(forward_match))
    print('Backward Match:', len(backward_match))

    df = df[df.index.isin(forward_match) | df.index.isin(backward_match)]

    if len(df) == 0:
        continue

    # Set the new variant ID according to the bim file
    df.loc[forward_match, 'SNP'] = bim1.loc[forward_match, 1]
    df.loc[backward_match, 'SNP'] = bim2.loc[backward_match, 1]
    ss_dfs.append(df)

try:
    df = pd.concat(ss_dfs).reset_index()
except ValueError:
    df = pd.DataFrame(columns=['CHR', 'BP', 'SNP', 'P'])

print(df)

# If p2 thresh is 1, then we want to keep all available SNPs
# This means that there could be variants in the ref panel which aren't in the GWAS
# So, we add new rows to the GWAS summary stats for those missing variants
# This comes in handy so we aren't dropping proxy SNPs that may be associated with other phenotypes

if p2thresh == 1 and len(df) > 0 and df['P'].min() < p1thresh:

    bim_non_match = bim[~bim[1].isin(df['SNP'])]
    bim_non_match[3] = bim_non_match[3].astype(int)

    # Pull out min/max windows
    candidate_lead_snp_pos = df[df['P'] <= p1thresh]['BP'].astype(int).sort_values()
    candidate_lead_snp_pos_min = candidate_lead_snp_pos - bp_window
    candidate_lead_snp_pos_max = candidate_lead_snp_pos + bp_window

    overlap_test = candidate_lead_snp_pos_min.iloc[1:].values > candidate_lead_snp_pos_max.iloc[:-1].values
    print(overlap_test, len(candidate_lead_snp_pos_min))

    while np.any(overlap_test):

        merge_row1 = candidate_lead_snp_pos_max.index[:-1][overlap_test][0]
        merge_row2 = candidate_lead_snp_pos_min.index[1:][overlap_test][0]

        new_start = min(candidate_lead_snp_pos_min.loc[merge_row1], candidate_lead_snp_pos_min.loc[merge_row2])
        new_end = max(candidate_lead_snp_pos_max.loc[merge_row1], candidate_lead_snp_pos_max.loc[merge_row2])

        candidate_lead_snp_pos_min = candidate_lead_snp_pos_min.drop(merge_row2)
        candidate_lead_snp_pos_max = candidate_lead_snp_pos_max.drop(merge_row2)

        candidate_lead_snp_pos_min.loc[merge_row1] = new_start
        candidate_lead_snp_pos_max.loc[merge_row1] = new_end

        overlap_test = candidate_lead_snp_pos_min.iloc[1:].values > candidate_lead_snp_pos_max.iloc[:-1].values
        print(overlap_test, len(candidate_lead_snp_pos_min))


    print(len(candidate_lead_snp_pos))

    bim_non_match = bim_non_match[(candidate_lead_snp_pos_min.min() < bim_non_match[3]) & (bim_non_match[3] < candidate_lead_snp_pos_max.max())]

    bim_in_window_mask = np.zeros(len(bim_non_match)).astype(bool)

    for start, end in zip(candidate_lead_snp_pos_min.values, candidate_lead_snp_pos_max.values):
        bim_in_window_mask |= ((start < bim_non_match[3]) & (bim_non_match[3] < end))

    add_rows = bim_non_match[bim_in_window_mask].copy()

    add_df = add_rows.rename(columns={0: 'CHR', 3: 'BP', 1: 'SNP'})[['CHR', 'BP', 'SNP']]
    add_df['P'] = 0.999

    df = pd.concat([df, add_df])


df = df.drop_duplicates(subset='SNP', keep='first')
df[['SNP', 'P', 'CHR', 'BP']].to_csv(output_sumstats, sep='\t', index=False, na_rep='NA')

# This is for the extract file
df['LABEL'] = df['SNP']
df[['CHR', 'BP', 'BP', 'LABEL']].dropna(how='any').to_csv(output_extract, sep='\t', index=False, header=False)

# This gets the min p value so we know whether to actually test this chromosome
p_series = pd.Series(dtype=object, index=['ANALYSIS', 'CHR', 'MIN_P'])
p_series.loc['ANALYSIS'] = analysis
p_series.loc['CHR'] = chr
p_series.loc['MIN_P'] = df['P'].min() if len(df) > 0 else 2
p_series.to_csv(output_min_p, header=False)