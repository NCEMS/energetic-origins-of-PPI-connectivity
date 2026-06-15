### `18-finalize`: Finally, produce the analysis-ready data product

Corresponding configuration file section: `FINALIZE NODES`

Anticipated execution time: ~1 min

To run this pipeline step in isolation, run the command:

```bash
snakemake -c all --use-conda --conda-frontend conda --snakefile 18-finalize/Snakefile --configfile config-files/s288c.config
```

* Write a few different "final" versions of the annotated protein-protein interation network
* As its name implies, this pipeline step produces a version of the database which is "flattened" over IDRs, with one row per IDR rather than one row per protein
