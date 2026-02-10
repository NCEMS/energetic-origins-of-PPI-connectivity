### `16-domain-annotations`: Integrate domain annotations from InterPro

Corresponding configuration file section: `Integrate domain annotations`

Anticipated execution time: <10 min

To run this pipeline step in isolation, run the command:

```bash
snakemake -c all --use-conda --conda-frontend conda --snakefile 14-domain-annotations/Snakefile --configfile config-files/s288c.config
```

* Parses and integrates domain annotations from the InterPro database
* Parsing the database file will take ~10-20 minutes on a typical computer
