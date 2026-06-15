# writes an empty PTM CSV with a header so downstream merges don't break.
# access params via 'snakemake' object provided by Snakemake.
from pathlib import Path
import pandas as pd

out = Path(str(snakemake.output[0]))
header = snakemake.params.get("header", "Accession,PTM_type,Description,Begin,End")
cols = [c.strip() for c in header.split(",") if c.strip()]

out.parent.mkdir(parents=True, exist_ok=True)
pd.DataFrame(columns=cols).to_csv(out, index=False)
