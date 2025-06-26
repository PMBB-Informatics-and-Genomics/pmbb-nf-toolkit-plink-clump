import pandas as pd
import sys

gwas = pd.read_table(sys.argv[1])
print(len(gwas), 'rows to start')
gwas = gwas.dropna(subset=['mapped_trait', 'mapped_uri'], how='any')
print(len(gwas), 'rows after dropping missing trait and URI values')
gwas = gwas[gwas['P'] <= 5E-8]
print(len(gwas), 'rows after dropping sub-significant P-values')

"""
gwas['mapped_trait'] = gwas['mapped_trait'].str.replace('osteoarthritis, ', 'osteoarthritis of ')
gwas['mapped_trait'] = gwas['mapped_trait'].str.replace('migraine without aura, susceptibility to, 4', 'inherited susceptibility to migraine without aura')
gwas['mapped_trait'] = gwas['mapped_trait'].str.replace('polyarticular juvenile idiopathic arthritis, rheumatoid factor negative', 'rheumatoid factor negative polyarticular juvenile idiopathic arthritis')
"""

gwas['num_traits'] = gwas['mapped_trait'].str.count(',') + 1
gwas['num_URIs'] = gwas['mapped_uri'].str.count(',') + 1

mismatch_fix = gwas.index[(gwas['num_traits'] != gwas['num_URIs']) & (gwas['num_URIs'] == 1)]
gwas.loc[mismatch_fix, 'mapped_trait'] = gwas.loc[mismatch_fix, 'mapped_trait'].str.replace(', ', ' - ')

gwas['num_traits'] = gwas['mapped_trait'].str.count(',') + 1
print(gwas[gwas['num_traits'] > 1])

new_rows = [gwas[gwas['num_traits'] == 1].transpose()]

for _, row in gwas[gwas['num_traits'] > 1].iterrows():
    if row['num_traits'] == 1:
        new_rows.append(row)
        continue
    trait_list = [t.strip() for t in row['mapped_trait'].split(',')]
    uri_list = [u.strip() for u in row['mapped_uri'].split(',')]

    subDF = pd.concat([row] * row['num_traits'], axis=1)
    subDF.loc['mapped_trait'] = trait_list
    subDF.loc['mapped_uri'] = uri_list
    new_rows.append(subDF)

new_df = pd.concat(new_rows, axis=1).transpose()
new_df = new_df.drop(columns=['num_traits', 'num_URIs'])
print(new_df)

new_df.to_csv(sys.argv[1] + "_processed.bed", sep='\t', index=False)