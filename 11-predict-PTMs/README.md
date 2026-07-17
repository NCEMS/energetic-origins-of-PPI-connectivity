### `10-predict-PTMs`: Predict post-translation modifications

Corresponding configuration file section: `Predict post-translational modifications`

```bash
micromamba env create -f env/ptmgpt2-v2.yml -n test-PTMGPT2-gpu --override-channels -y
micromamba activate test-PTMGPT2-gpu
micromamba --override-channels -c pytorch -c nvidia -c conda-forge -c bioconda snakemake -y
snakemake -c 96 --configfile config-files/union-athaliana.config --snakefile 11-predict-PTMs/Snakefile
```

### Legacy yeast information follows below


Anticipated execution time: 2-3 days (full run with 2 x GPUs) or 3 min (reuse pre-computed PTMs)

To run this pipeline step in isolation, run the command:

```bash
snakemake -c all --use-conda --conda-frontend conda --snakefile 10-predict-PTMs/Snakefile --configfile config-files/s288c.config
```

If you do not want to rerun predictions, use the command:

```bash
snakemake -c 112 --use-conda --conda-frontend conda --configfile config-files/s288c.config --snakefile 10-predict-PTMs/Snakefile --allowed-rules process_PTMs
```

This command assumes that you have downloaded the pre-computed results from CyVerse and unpacked them into processed-data as described in `../README.md`

* This step runs the model PTMGPT2 to predict the locations of 19 different post-translational modifications within each protein sequence
* Predictions are run on SignalP-trimmed sequences, but residue indices are mapped back to the original numbering for the untruncated sequence. 
* The current Snakefile assumes you have 2 GPUs with device IDs [0, 1]. Using two RTX6000 GPUs in coarse-grain parallel, predictions required ~48 hours for all 19 models. If you have a single GPU only, you can run all predictions on one GPU by replacing `run-PTMGPT2-parallel.py` on line 39 of the `Snakefile` (within the rule `predict_PTMs`) with `run-PTMGPT2.py`, which will look for a single GPU only. 
