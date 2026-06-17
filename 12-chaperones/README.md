### `12-chaperones`: Integrate information on chaperone interactions

```text
Input nodes: 11,695
Unique chaperone loci in chaperone table: 305
Nodes with at least one interacting chaperone: 1,122
Top chaperone family counts among annotated node interactions:
interacting_chaperone_families
CYP-type PPIases               448
Chaperonins                    189
HSP90 co-chaperones            142
Small HSPs                     125
HSP90s                         117
HSP70s                         112
DNAJ proteins                   97
FKBP-type PPIases               83
Protein disulfide isomerase     65
HSP100s                         18
HSP110s                         14
BAGs                            11
```

### LEGACY YEAST INFORMATION FOLLOWS

Corresponding configuration file section: `Integrate chaperone information`

Anticipated execution time: 2 min

To run this pipeline step in isolation, run the command:

```bash
snakemake -c all --use-conda --conda-frontend conda --snakefile 12-chaperones/Snakefile --configfile config-files/s288c.config
```

* This pipeline section uses a list of unique identifiers for chaperones and co-chaperones in yeast to determine which proteins in The Yeast Interactome interact with which chaperones
* The network edges from your input interaction network in `1-networt-centrality` will be used to determine chaperone interactions
