### `9-translation-efficiency`: Incorporate translation efficiency data

Corresponding configufle section: `Translation efficiency`

Anticipated execution time: 1 min

To run this pipeline step in isolation, run the command:

```bash
snakemake -c all --use-conda --conda-frontend conda --snakefile 9-translation-efficiency/Snakefile --configfile config-files/s288c.config
```


* Merges translation efficiency information from `scikit-ribo` into the dataset
* Original data are contained within the `scikit-ribo_manuscript` repo [here](https://github.com/schatzlab/scikit-ribo_manuscript/blob/master/Data/skr_weinberg_genesTE.csv). 
