### `7-protein-half-life`: Incorporate protein half life data

Corresponding configfile section: `Protein half-life`

Anticipated execution time: 1 min

In the Fan2016 dataset (DB1), the column for half-life is labelled with units of "hr-1", which I am assuming is actually hours. 

### Legacy s288c information

To run this pipeline step in isolation, run the command:

```bash
snakemake -c all --use-conda --conda-frontend conda --snakefile 7-protein-half-life/Snakefile --configfile config-files/s288c.config
```

* Merges protein half-life data into the dataset
* Data from [Christiano et al. 2014](https://doi.org/10.1016/j.celrep.2014.10.065) and [Martin-Perez and Villen](https://doi.org/10.1016/j.cels.2017.08.008) are both integrated
