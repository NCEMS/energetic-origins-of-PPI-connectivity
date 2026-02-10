### `8-protein-expression`: Incorporate protein expression data

Corresponding configfile section: `Protein expression`

Anticipated execution time: 1 min

To run this pipeline step in isolation, run the command:

```bash
snakemake -c all --use-conda --conda-frontend conda --snakefile 8-protein-expression/Snakefile --configfile config-files/s288c.config
```


* Merges protein expression data into the dataset
* The default .config file will merge in data from `../0-download-inputs/data-files/1-s2.0-S240547121730546X-mmc5.csv`, which is from Table S4 of 10.1016/j.cels.2017.12.004

Note well: this pipeline depends on specific column names existing within the expression data to be merged into the protein interaction network. If you want to use different data and have a new file, consult the contents of `python-scripts/add-expression.py` and consider if the column names used for merging need to be adjusted. 
