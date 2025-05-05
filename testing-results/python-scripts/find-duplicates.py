import os, sys
import pandas as pd

# load a nodes_df .csv file
df = pd.read_csv(sys.argv[1])

# group and filter to find duplicates
filtered = df.groupby('UniProtKB-AC').filter(lambda g: g['node'].nunique() > 1)

# print summary information
print (filtered.info())
print (filtered.head())

# extract the two columns required
temp = filtered[["UniProtKB-AC", "node"]].copy()

# nodes
print (temp)

# UniProt IDs
#print (filtered["UniProtKB-AC"])
