#!/usr/bin/env python3

import pandas as pd
import pint
import pint_pandas

file_path = "18-entanglement/processed-data/union-athaliana-step18.pkl"

df = pd.read_pickle(file_path)

print(f"Number of columns: {len(df.columns)}\n")

for column in df.columns:
    print(column)
