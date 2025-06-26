#!/usr/bin/env python3

import pandas as pd
import argparse

def main():
    parser = argparse.ArgumentParser(description='Create fake summary statistics')
    parser.add_argument('--bim', required=True, help='BIM file path')
    parser.add_argument('--snps', required=True, help='Lead SNp_values file path')
    parser.add_argument('--output', required=True, help='Output file path')
    
    args = parser.parse_args()
    
    # Read input files
    snps = pd.read_csv(args.snps)
    print(snps)
    snps = snps.set_index(['chromosome', 'base_pair_location'])
    print(snps)
    
    bim = pd.read_table(args.bim, header=None, names=['chromosome', 'variant_id', 'CM', 'base_pair_location', 'other_allele', 'effect_allele'])
    bim = bim.set_index(['chromosome', 'base_pair_location'])
    print(bim[bim.index.isin(snps.index)])
    
    # Assign p_value-values
    bim['p_value'] = 0.9
    bim.loc[snps.index.intersection(bim.index), 'p_value'] = 0.001
    bim = bim[['variant_id', 'p_value', 'other_allele', 'effect_allele']].reset_index()
    print(bim)
    
    # Write output
    bim.to_csv(args.output, sep='\t', index=False)

if __name__ == "__main__":
    main()