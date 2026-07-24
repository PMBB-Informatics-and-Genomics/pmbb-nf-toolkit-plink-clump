import pandas as pd
import numpy as np
import json
import sys

#!/usr/bin/env python3
import sys
import argparse

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser()
    
    parser.add_argument('--clumps', required=True, help='Path to the clump summary file')
    parser.add_argument('--loci', required=True, help='Path to the locus summary file')
    parser.add_argument('--gwas', required=True, help='Path to the GWAS catalog file')
    parser.add_argument('--traits', required=False, nargs='*')
    
    return parser.parse_args()

args = parse_arguments()
    
clumpFile = args.clumps
locusFile = args.loci
gwasCatalogFile = args.gwas
traits_of_interest = args.traits if args.traits is not None else []

clumps = pd.read_csv(clumpFile)
clumps[['#CHROM', 'POS']] = clumps[['#CHROM', 'POS']].astype(int)
clumps = clumps.rename(columns={'#CHROM': 'CHR', 'POS': 'BP'})
print(clumps) 

loci = pd.read_csv(locusFile)
loci[['#CHROM', 'POS']] = loci[['#CHROM', 'POS']].astype(int)
loci = loci.rename(columns={'#CHROM': 'CHR', 'POS': 'BP'})
print(loci) 

gwas_all = pd.read_table(gwasCatalogFile)
gwas_all['CHR'] = pd.to_numeric(gwas_all['CHR_ID'].str.replace('chr', '').replace({'X': 23, 'Y': 24}), errors='coerce')
gwas_all['BP'] = pd.to_numeric(gwas_all['CHR_POS'], errors='coerce')
gwas_all = gwas_all.dropna(subset=['CHR', 'BP'])
gwas_all = gwas_all.set_index(['CHR', 'BP'], drop=False).sort_index()
print(gwas_all)
print(gwas_all.columns)

for result_type, results_df in {'clumps': clumps, 'loci': loci}.items():
    for idx, row in results_df.iterrows():
        coords = (row['CHR'], row['BP'])
        min_pos, max_pos = row['MIN_POS'], row['MAX_POS']

        if coords in gwas_all.index:
            direct_match_rows = gwas_all.loc[coords]
            results_df.loc[idx, 'Lead_SNP_Match_With_Previosuly_Known'] = True
            matches = [
                (
                    g_row['MAPPED_TRAIT'],
                    g_row['MAPPED_TRAIT_URI'].split('/')[-1],
                    f'PMID:{g_row["PUBMEDID"]}',
                    g_row['STUDY ACCESSION']
                ) for _, g_row in direct_match_rows.iterrows()
            ]
            matches = list(set(matches))
            results_df.loc[idx, 'Lead_SNP_Matches_JSON'] = json.dumps(matches)
        else:
            results_df.loc[idx, 'Lead_SNP_Match_With_Previosuly_Known'] = False
        
        gwas_match = gwas_all[(gwas_all['CHR'] == row['CHR']) & (gwas_all['BP'].between(min_pos, max_pos))]
        gwas_match = gwas_match[gwas_match.index != coords]

        if len(gwas_match) == 0:
            results_df.loc[idx, 'Region_Overlap_With_Previously_Known'] = False | results_df.loc[idx, 'Lead_SNP_Match_With_Previosuly_Known']
            continue
        else:
            results_df.loc[idx, 'Region_Overlap_With_Previously_Known'] = True
            matches = [
                (
                    g_row['MAPPED_TRAIT'],
                    g_row['MAPPED_TRAIT_URI'].split('/')[-1],
                    f'{int(g_row["CHR"])}:{int(g_row["BP"])}',
                    f'PMID:{g_row["PUBMEDID"]}',
                    g_row['STUDY ACCESSION']
                ) for _, g_row in gwas_match.iterrows()
            ]
            matches = list(set(matches))
            results_df.loc[idx, 'Region_Overlap_Matches_JSON'] = json.dumps(matches)
    
    results_df[['Lead_SNP_Match_With_Previosuly_Known', 'Region_Overlap_With_Previously_Known']] = results_df[['Lead_SNP_Match_With_Previosuly_Known', 'Region_Overlap_With_Previously_Known']].astype(bool)

    if len(traits_of_interest) > 0:
        results_df['Lead_SNP_Match_In_Traits_Of_Interest'] = results_df['Lead_SNP_Matches_JSON'].astype(str).apply(lambda x: np.any([t.lower() in x.lower() for t in traits_of_interest]))
        results_df['Region_Overlap_Match_In_Traits_Of_Interest'] = results_df['Region_Overlap_Matches_JSON'].astype(str).apply(lambda x: np.any([t.lower() in x.lower() for t in traits_of_interest]))
        results_df['Region_Overlap_Match_In_Traits_Of_Interest'] |= results_df['Lead_SNP_Match_In_Traits_Of_Interest']

        results_df[['Lead_SNP_Match_In_Traits_Of_Interest', 'Region_Overlap_Match_In_Traits_Of_Interest']] = results_df[['Lead_SNP_Match_In_Traits_Of_Interest', 'Region_Overlap_Match_In_Traits_Of_Interest']].astype(bool)

    print(results_df)

    with open(f'gwas_catalog_{result_type}_lookup_count_summary.txt', 'w+') as f:
        print(f'Total number of {result_type}:', len(results_df), '\n', file=f)
        print('---Summary counts:', file=f)
        if len(traits_of_interest) > 0:
            print(results_df[[
                'Region_Overlap_With_Previously_Known',
                'Lead_SNP_Match_With_Previosuly_Known',
                'Region_Overlap_Match_In_Traits_Of_Interest',
                'Lead_SNP_Match_In_Traits_Of_Interest'
            ]].apply(lambda x: x.value_counts(), result_type='expand').transpose(), '\n', file=f)
        else:
            print(results_df[[
                'Region_Overlap_With_Previously_Known',
                'Lead_SNP_Match_With_Previosuly_Known'
            ]].apply(lambda x: x.value_counts(), result_type='expand').transpose(), '\n', file=f)

        print('---Any trait results:', file=f)
        print(f'{result_type.capitalize()} with no matches in the GWAS Catalog:', (~results_df['Region_Overlap_With_Previously_Known']).sum(), file=f)
        print(f'{result_type.capitalize()} with any matches in the GWAS Catalog:', (results_df['Region_Overlap_With_Previously_Known']).sum(), file=f)
        print(f'Of those, {result_type} with lead SNP matches in the GWAS Catalog:', (results_df['Lead_SNP_Match_With_Previosuly_Known']).sum(), '\n', file=f)

        if len(traits_of_interest) > 0:
            print('---Traits of interest results:', file=f)
            print(f'{result_type.capitalize()} with no trait of interest matches in the GWAS Catalog:', (~results_df['Region_Overlap_Match_In_Traits_Of_Interest']).sum(), file=f)
            print(f'{result_type.capitalize()} with trait of interest matches in the GWAS Catalog:', (results_df['Region_Overlap_Match_In_Traits_Of_Interest']).sum(), file=f)
            print(f'Of those, {result_type} with lead SNP matches with traits of interest in the GWAS Catalog:', (results_df['Lead_SNP_Match_In_Traits_Of_Interest']).sum(), '\n', file=f)

    results_df.to_csv(f'all_{result_type}.gwas_catalog_overlap.csv', index=False)

