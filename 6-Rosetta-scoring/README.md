### `6-Rosetta-scoring`: Run Rosetta FastRelax and scoring

Corresponding configuration file section: `Rosetta scoring`

* See `../README.md` for information on how to download the pre-relaxed structures and score files from CyVerse
* Uses the protein structures from AF2 as inputs to a Rosetta FastRelax protocol with the ref2015 scoring function
* You must update your .config file to include the absolute path the your Rosetta `relax.static.linuxgccrelease` executable in the corresponding section
* In the current implementation, 10 independent replicas were run for each protein
* Note that this step take significant time; uses Python multiprocessing to distribute calculations over N CPUs as specified in the .config file for the run.
* With N = 96 CPUs, this step required ~30 days. 
