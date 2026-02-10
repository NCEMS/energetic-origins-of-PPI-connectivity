### `19-meltome-atlas`: Integrate protein stability information from the Meltome Atlas

Corresponding configuration file section: `Integrate TPP data from Meltome Atlas`

Anticipated execution time: ~5 min

To run this pipeline step in isolation, run the command:

```bash
snakemake -c all --use-conda --conda-frontend conda --snakefile 17-meltome-atlas/Snakefile --configfile config-files/s288c.config
```

* Meltome Atlas data were extracted from Table S2 of [10.1038/s41592-020-0801-4](https://www.nature.com/articles/s41592-020-0801-4), specifically the tab labelled `ma_0010`, corresponding to the sample `Saccharomyces_cerevisiae_20171013`
