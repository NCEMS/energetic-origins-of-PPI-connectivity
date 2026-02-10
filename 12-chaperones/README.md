### `14-chaperones`: Integrate information on chaperone interactions

Corresponding configuration file section: `Integrate chaperone information`

Anticipated execution time: 2 min

To run this pipeline step in isolation, run the command:

```bash
snakemake -c all --use-conda --conda-frontend conda --snakefile 12-chaperones/Snakefile --configfile config-files/s288c.config
```

* This pipeline section uses a list of unique identifiers for chaperones and co-chaperones in yeast to determine which proteins in The Yeast Interactome interact with which chaperones
* The network edges from your input interaction network in `1-networt-centrality` will be used to determine chaperone interactions
