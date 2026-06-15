### `5-dG-calculations`: Predict protein stability by a few methods

Corresponding configuration file section: `ΔG PREDICTIONS`

Anticipated execution time: <30 min (with CUDA)

To run this pipeline step in isolation, run the command:

```bash
snakemake -c all --use-conda --conda-frontend conda --snakefile 5-dG-calculations/Snakefile --configfile config-files/s288c.config
```

* Prepare protein structures for dG predictions by locating them in the pipeline directories and creating cleaved PDBs as needed (logging information in `prepare-structures.log`)
* Use Eq. 1 from Ghosh & Dill 2010 to predict dG
* Use Cagiada ESM-IF based LLM to predict dG (logging information in `cagiada.log`)

* N.B. - this step can be directed towards a folder with new AlphaFold2 predictions to supplement those from EBI. If you do not have new predictions to fill in structures with a cleavage site or sequence mismatch, you can direct the pipeline to a dummy directory. No "rescue" structures will be found, so the standard EBI structures will be used as possible. 
* If you want to use the new structure predictions, download the required data from CyVerse:

```bash
gocmd get --progress /iplant/home/shared/NCEMS/working-groups/energetic-origins/required-data/5-dG-calculations/alphafold2_structures /your/selected/path/for/structures
```

* Update the config variable "AF2_dir" to point to `{your-local-path}/alphafold2-structures`
