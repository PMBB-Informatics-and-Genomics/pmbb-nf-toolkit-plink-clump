import pandas as pd
import numpy as np
import json
import sys
from fuzzywuzzy import fuzz

def fuzzy_match(string, search_description):
    if not isinstance(string, str):
        string = str(string)
    return fuzz.token_sort_ratio(string.lower(), search_description.lower())

#!/usr/bin/env python3
import sys
import argparse

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='Process genomic data files.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument('--clump', required=True,
                        help='Path to the clump file')
    
    parser.add_argument('--snp', required=True,
                        help='Path to the SNP file')
    
    parser.add_argument('--gwas-catalog', required=True,
                        help='Path to the GWAS catalog file')
    
    parser.add_argument('--trait-map', required=True,
                        help='Path to the trait map file')
    
    parser.add_argument('--chunkmap', required=True,
                        help='Path to the chunkmap file')
    
    parser.add_argument('--coords', required=True,
                        help='Path to the coordinates file')
    
    parser.add_argument('--trait-map-out', required=True,
                        help='Path to write the trait map output')
    
    parser.add_argument('--tab-out', required=True,
                        help='Path to write the tab output')
    
    return parser.parse_args()


args = parse_arguments()
    
clumpFile = args.clump
snpFile = args.snp
gwasCatalogFile = args.gwas_catalog
traitMapFile = args.trait_map
chunkmapFile = args.chunkmap
coordsFile = args.coords
traitmapOut = args.trait_map_out
tabOut = args.tab_out

clumps = pd.read_table(clumpFile, sep='\s+')
clumps[['CHR', 'BP']] = clumps[['CHR', 'BP']].astype(int)
clumps = clumps.set_index(['CHR', 'BP'])
print(clumps) 

gwas = pd.read_table(gwasCatalogFile)
gwas['CHR'] = gwas['chr'].str.replace('chr', '').replace({'X': 23, 'Y': 24}).astype(int)
gwas['BP'] = gwas['start']
gwas[['CHR', 'BP']] = gwas[['CHR', 'BP']].astype(int)
gwas = gwas.set_index(['CHR', 'BP']).sort_index()
print(gwas)

gwas_all = gwas.copy()

snps = pd.read_table(snpFile)
snps = snps.sort_values(by=['CHR', 'BP'])
snps[['CHR', 'BP']] = snps[['CHR', 'BP']].astype(int)
snps = snps.set_index(['CHR', 'BP'])
print(snps)

#think about casing later

print(gwas)

#terms provided by the user params
#search_terms = [u.split('/')[-1] for u in search_URI]
gwas_all = gwas_all.sort_index()
print(gwas_all)
gwas['trait_efo'] = gwas['MAPPED_TRAIT_URI'].str.split('/').apply(lambda x: x[-1])
trait_info = gwas[['trait_efo', 'MAPPED_TRAIT', 'MAPPED_TRAIT_URI']].drop_duplicates()

print(trait_info)
trait_info.to_csv(traitmapOut, sep='\t', index=False)

print(len(snps), len(clumps))
print(len(snps.index.intersection(clumps.index)))

clumps = clumps.reset_index()
snps['Novel'] = np.nan
snps['Novel_Logic'] = 'SNP is Not in Clumps'
snps['Novel_SNP'] = np.nan
snps['Novel_SNP_Logic'] = 'SNP is Not in Clumps'
snps['Novel_SNP_Coords'] = ~snps.index.isin(gwas_all.index)

direct_matches = snps.index.intersection(gwas.index)

snps.loc[direct_matches, 'Novel'] = False
snps.loc[direct_matches, 'Novel_Logic'] = 'Lead SNP found in GWAS Catalog'

direct_matches = snps.index.intersection(gwas_all.index)

snps.loc[direct_matches, 'Novel_SNP'] = False
snps.loc[direct_matches, 'Novel_SNP_Logic'] = 'Lead SNP found in GWAS Catalog'

for _, row in clumps.iterrows():
    coords = (row['CHR'], row['BP'])

    if row['SP2'] == 'NONE' and (coords not in gwas.index): # Lead SNP has no tagging SNPs
        snps.loc[coords, 'Novel'] = True
        snps.loc[coords, 'Novel_Logic'] = 'Lead SNP not in GWAS Catalog and Has no tagging SNPs'
        continue
    elif row['SP2'] == 'NONE' and (coords in gwas.index):
        snps.loc[coords, 'Novel'] = False
        snps.loc[coords, 'Novel_Logic'] = 'Lead SNP found in GWAS Catalog'
        continue

    # Extract tag SNPs from clump row
    tags = pd.Series(row['SP2'].split(','))
    tags = tags.str.split('(', expand=True)[0]
    tags = tags.str.split(':', expand=True)
    tags = tags.reset_index()
    tags[[0, 1]] = tags[[0, 1]].astype(int)
    tags = tags.set_index([0, 1]).sort_index()
    # tags[[0, 1]] = tags[[0, 1]].astype(int)
    # tags = tags.set_index([0, 1])

    if coords in gwas.index: # Lead SNP found in GWAS catalog
        snps.loc[coords, 'Novel'] = False
        snps.loc[coords, 'Novel_Logic'] = 'Lead SNP found in GWAS Catalog'
        # Also update the Novel value for tag SNPs if the lead SNP is in the GWAS catalog
        if len(snps.index.intersection(tags.index)) > 0:
            lead_tag_coords = snps.index.intersection(tags.index)
            snps.loc[lead_tag_coords, 'Novel'] = False
            snps.loc[lead_tag_coords, 'Novel_Logic'] = 'Lead SNP became Tag SNP, Lead SNP of Clump found in GWAS Catalog'
        continue

    # Test if any of the tag SNPs are in the GWAS catalog
    if len(tags.index.intersection(gwas.index)) == 0:
        # Lead SNP not found and tag SNPs not found
        snps.loc[coords, 'Novel'] = True
        snps.loc[coords, 'Novel_Logic'] = 'No Tag SNPs are in the GWAS Catalog'
        if len(snps.index.intersection(tags.index)) > 0:
            # Set all tag SNPs to novel if lead SNP and tags not found
            lead_tag_coords = snps.index.intersection(tags.index)
            snps.loc[lead_tag_coords, 'Novel'] = True
            snps.loc[lead_tag_coords, 'Novel_Logic'] = 'Lead SNP became Tag SNP, Other SNPs in Clump not found in GWAS Catalog'
    else:
        # Lead SNP not found but tag SNP(s) found
        snps.loc[coords, 'Novel'] = False
        snps.loc[coords, 'Novel_Logic'] = 'Tag SNP(s) found in GWAS Catalog'
        if len(snps.index.intersection(tags.index)) > 0:
            # Set all tag SNPs to known if lead SNP or tags found
            lead_tag_coords = snps.index.intersection(tags.index)
            snps.loc[lead_tag_coords, 'Novel'] = False
            snps.loc[lead_tag_coords, 'Novel_Logic'] = 'Lead SNP became Tag SNP, Other SNPs in Clump found in GWAS Catalog'

for _, row in clumps.iterrows():
    coords = (row['CHR'], row['BP'])

    if row['SP2'] == 'NONE' and (coords not in gwas_all.index): # Lead SNP has no tagging SNPs
        snps.loc[coords, 'Novel_SNP'] = True
        snps.loc[coords, 'Novel_SNP_Logic'] = 'Lead SNP not in GWAS Catalog and Has no tagging SNPs'
        continue
    elif row['SP2'] == 'NONE' and (coords in gwas_all.index):
        snps.loc[coords, 'Novel_SNP'] = False
        snps.loc[coords, 'Novel_SNP_Logic'] = 'Lead SNP found in GWAS Catalog'
        continue

    # Extract tag SNPs from clump row
    tags = pd.Series(row['SP2'].split(','))
    tags = tags.str.split('(', expand=True)[0]
    tags = tags.str.split(':', expand=True)
    tags = tags.reset_index()
    tags[[0, 1]] = tags[[0, 1]].astype(int)
    tags = tags.set_index([0, 1]).sort_index()
    # tags[[0, 1]] = tags[[0, 1]].astype(int)
    # tags = tags.set_index([0, 1])

    if coords in gwas_all.index: # Lead SNP found in GWAS catalog
        snps.loc[coords, 'Novel_SNP'] = False
        snps.loc[coords, 'Novel_SNP_Logic'] = 'Lead SNP found in GWAS Catalog'
        # Also update the Novel value for tag SNPs if the lead SNP is in the GWAS catalog
        if len(snps.index.intersection(tags.index)) > 0:
            lead_tag_coords = snps.index.intersection(tags.index)
            snps.loc[lead_tag_coords, 'Novel_SNP'] = False
            snps.loc[lead_tag_coords, 'Novel_SNP_Logic'] = 'Lead SNP became Tag SNP, Lead SNP of Clump found in GWAS Catalog'
        continue

    # Test if any of the tag SNPs are in the GWAS catalog
    if len(tags.index.intersection(gwas_all.index)) == 0:
        # Lead SNP not found and tag SNPs not found
        snps.loc[coords, 'Novel_SNP'] = True
        snps.loc[coords, 'Novel_SNP_Logic'] = 'No Tag SNPs are in the GWAS Catalog'
        if len(snps.index.intersection(tags.index)) > 0:
            # Set all tag SNPs to novel if lead SNP and tags not found
            lead_tag_coords = snps.index.intersection(tags.index)
            snps.loc[lead_tag_coords, 'Novel_SNP'] = True
            snps.loc[lead_tag_coords, 'Novel_SNP_Logic'] = 'Lead SNP became Tag SNP, Other SNPs in Clump not found in GWAS Catalog'
    else:
        # Lead SNP not found but tag SNP(s) found
        snps.loc[coords, 'Novel_SNP'] = False
        snps.loc[coords, 'Novel_SNP_Logic'] = 'Tag SNP(s) found in GWAS Catalog'
        if len(snps.index.intersection(tags.index)) > 0:
            # Set all tag SNPs to known if lead SNP or tags found
            lead_tag_coords = snps.index.intersection(tags.index)
            snps.loc[lead_tag_coords, 'Novel_SNP'] = False
            snps.loc[lead_tag_coords, 'Novel_SNP_Logic'] = 'Lead SNP became Tag SNP, Other SNPs in Clump found in GWAS Catalog'

#is known association is with phenotype terms of interest or else
#cite the study accession in output table
print(snps['Novel'].fillna('NA').value_counts(), '\n')
print(snps['Novel_Logic'].value_counts(), '\n')
print(snps['Novel_SNP'].fillna('NA').value_counts(),'\n')
print(snps['Novel_SNP_Logic'].value_counts(),'\n')
print(snps['Novel_SNP_Coords'].value_counts(), '\n')

print('\nNew Columns\n')
snps['Known_Association'] = np.logical_not(snps['Novel'])
snps['Known_Association'] = snps['Known_Association'].mask(pd.isnull(snps['Novel']))

snps['Novel_Association_Known_Signal'] = snps['Novel'] & np.logical_not(snps['Novel_SNP'])
snps['Novel_Association_Known_Signal'] = snps['Novel_Association_Known_Signal'].mask(pd.isnull(snps['Novel']) | pd.isnull(snps['Novel_SNP']))

snps['Novel_Signal'] = snps['Novel_SNP']

print(snps['Known_Association'].fillna('NA').value_counts(),'\n')
print(snps['Novel_Association_Known_Signal'].fillna('NA').value_counts(),'\n')
print(snps['Novel_Signal'].fillna('NA').value_counts(),'\n')

print(snps[['Known_Association', 'Novel_Association_Known_Signal', 'Novel_Signal']].fillna('NA').value_counts(),'\n')

snps.to_csv(tabOut+'.tsv', sep='\t', na_rep='NA') 
