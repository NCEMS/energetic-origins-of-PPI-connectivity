### `11-entanglement`: Integrate entanglement information into the interactome

Corresponding configuration file section: `Integrate entanglements from O'Brien Lab`

Anticipated execution time: 1 min

To run this pipeline step in isolation, run the command:

```bash
snakemake -c all --use-conda --conda-frontend conda --snakefile 11-entanglement/Snakefile --configfile config-files/s288c.config
```

* This pipeline annotations proteins with whether or not their experimental and/or AlphaFold2 structures contain entanglements
