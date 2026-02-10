### `15-oligomers`: Integrate information about known oligomerization states for proteins

Corresponding configuration file section: `Integrate oligomer information`

Anticipated execution time: 1 min

To run this pipeline step in isolation, run the command:

```bash
snakemake -c all --use-conda --conda-frontend conda --snakefile 13-oligomers/Snakefile --configfile config-files/s288c.config
```

* Integrate information about protein complexes from the Complex Portal into the interaction network
* Proteins in complexes with both known and unknown stoichiometry are considered separately
