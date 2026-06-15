### `6-Rosetta-scoring`: Run Rosetta FastRelax and scoring

Corresponding configuration file section: `Rosetta scoring`

Anticipated execution time: 40 days (full run with 112 CPUs) or 2 min (to reuse Rosetta calculations)

To run this pipeline step in isolation, run the command:

```bash
snakemake -c all --use-conda --conda-frontend conda --snakefile 6-Rosetta-scoring/Snakefile --configfile config-files/s288c.config
```

To reuse Rosetta structures and scores from CyVerse, target the final rule only with the command:

```bash
touch 6-Rosetta-scoring/processed-data/scores/.all_scores_done
snakemake -c all --use-conda --conda-frontend conda --snakefile 6-Rosetta-scoring/Snakefile --configfile config-files/s288c.config -- add_scores_to_nodes
```

* See `../README.md` for information on how to download the pre-relaxed structures and score files from CyVerse
* Uses the protein structures from AF2 as inputs to a Rosetta FastRelax protocol with the ref2015 scoring function
* You must update your .config file to include the absolute path the your Rosetta `relax.static.linuxgccrelease` executable in the corresponding section
* In the current implementation, 10 independent replicas were run for each protein
* Note that this step take significant time; uses Python multiprocessing to distribute calculations over N CPUs as specified in the .config file for the run.
* With N = 96 CPUs, this step required ~30 days. 
