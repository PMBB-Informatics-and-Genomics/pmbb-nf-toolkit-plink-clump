import pandas as pd
import sys

gwas = pd.read_table(sys.argv[1])
print(gwas)
print(len(gwas), 'rows to start')
gwas = gwas.dropna(subset=['MAPPED_TRAIT', 'MAPPED_TRAIT_URI'], how='any')
print(len(gwas), 'rows after dropping missing trait and URI values')
gwas = gwas[gwas['P-VALUE'] <= 5E-8]
print(len(gwas), 'rows after dropping sub-significant P-values')

# Now we can cound the number of traits represented by each row
gwas['num_URIs'] = gwas['MAPPED_TRAIT_URI'].str.count(',') + 1

new_rows = [gwas[gwas['num_URIs'] == 1].transpose()]

for _, row in gwas[gwas['num_URIs'] > 1].iterrows():
    # Split URIs by commas
    uri_list = [u.strip() for u in row['MAPPED_TRAIT_URI'].split(',')]

    # Split mapped traits by commas with max split to match length.
    # Might break trait names with commans in them!
    trait_list = [t.strip() for t in row['MAPPED_TRAIT'].split(',', len(uri_list)-1)]

    subDF = pd.concat([row] * row['num_URIs'], axis=1)
    subDF.loc['MAPPED_TRAIT', :] = trait_list
    subDF.loc['MAPPED_TRAIT_URI', :] = uri_list
    new_rows.append(subDF)

new_df = pd.concat(new_rows, axis=1).transpose()
new_df = new_df.drop(columns=['num_URIs'])
print(new_df)

new_df.to_csv("gwas_catalog_processed_with_split_rows.tsv.gz", sep='\t', index=False)