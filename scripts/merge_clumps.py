import sys
import pandas as pd

input_files = sys.argv[1:-1]
output_file = sys.argv[-1]

dfs = []
for f in input_files:
    try:
        dfs.append(pd.read_table(f, sep='\s+'))
    except pd.errors.EmptyDataError:
        continue

if len(dfs) == 0:
    cols = ['CHR','F','SNP','BP','P','TOTAL','NSIG','S05','S01','S001','S0001','SP2']
    results = pd.DataFrame(columns=cols)
else:
    results = pd.concat(dfs)

results.to_csv(output_file, sep='\t', index=False)