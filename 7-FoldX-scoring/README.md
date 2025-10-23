### `7-FoldX-scoring`: Run FoldX scoring 

Corresponding configfile section: `FoldX scoring`

* See `../README.md` for information on how to download pre-computed FoldX scores from CyVerse
* Runs FoldX scoring on Rosetta-relaxed poses; symlinks will be used to create local PDBs on which to run FoldX
* Utilizes Python multiprocessing to run calculations in parallel; you can choose the number of processed used in the config file
