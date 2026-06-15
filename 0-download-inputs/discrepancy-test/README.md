# Arabidopsis PPI Curation Discrepancy Test

This directory is an isolated scratch area for checking whether
`curate-athaliana-ppi-evidence.py` regenerates the existing manually curated
Arabidopsis PPI evidence tables.

It does not write into `../data-files`.

## Inputs

Place the raw, uncommitted source files here:

```text
inputs/BIOGRID-ARABIDOPSIS-THALIANA-5.0.258.tab3.xlsx
inputs/Table S1.xlsx
```

The runner also accepts explicit paths, so the raw inputs can stay elsewhere.

## Run

From this directory:

```bash
python run-discrepancy-test.py
```

Or from `0-download-inputs`:

```bash
python discrepancy-test/run-discrepancy-test.py
```

If your raw BioGRID file is a `.txt`, `.tab`, `.tsv`, or `.csv` file, pass it
explicitly:

```bash
python run-discrepancy-test.py \
  --biogrid path/to/BIOGRID-ARABIDOPSIS-THALIANA-5.0.258.tab3.txt \
  --table-s1 "path/to/Table S1.xlsx"
```

If starting from the full BioGRID ALL TAB3 file, first create an
Arabidopsis-only slice:

```bash
python filter-biogrid-athaliana.py \
  "path/to/BIOGRID-ALL-5.0.258.tab3.txt" \
  inputs/BIOGRID-ARABIDOPSIS-THALIANA-5.0.258.tab3.txt
```

Then run the discrepancy test using the manual provenance label:

```bash
python run-discrepancy-test.py \
  --biogrid inputs/BIOGRID-ARABIDOPSIS-THALIANA-5.0.258.tab3.txt \
  --table-s1 "path/to/Table S1.xlsx" \
  -- \
  --biogrid-source-file-label BIOGRID-ARABIDOPSIS-THALIANA-5.0.258.tab3.xlsx
```

## Outputs

All generated files are written to `outputs/`:

```text
generated_ppi_in_vivo_edges.csv
generated_ppi_in_vitro_edges.csv
unknown_methods.csv
discrepancy_summary.txt
in_vivo_only_manual_exact.csv
in_vivo_only_generated_exact.csv
in_vivo_keys_only_manual.csv
in_vivo_keys_only_generated.csv
in_vivo_metadata_mismatch.csv
in_vitro_only_manual_exact.csv
in_vitro_only_generated_exact.csv
in_vitro_keys_only_manual.csv
in_vitro_keys_only_generated.csv
in_vitro_metadata_mismatch.csv
```

The most important distinction:

* `*_keys_only_*` means the interaction pair itself differs.
* `*_metadata_mismatch.csv` means the pair exists in both files, but evidence
  count/source/method text differs.
