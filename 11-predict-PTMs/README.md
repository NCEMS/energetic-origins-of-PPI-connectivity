### `11-predict-PTMs`: Predict post-translation modifications

* Runs the model PTMGPT2 to predict the locations of 19 different post-translational modifications within each protein sequence

* Predictions are run on SignalP-trimmed sequences, but residue indices are mapped back to the original numbering for the untruncated sequence. 

*Note*: The current Snakefile assumes you have 2 GPUs with device IDs [0, 1]. Using two RTX6000 GPUs in coarse-grain parallel, predictions required ~48 hours for all 19 models. If you have a single GPU only, you can run all predictions on one GPU by replacing `run-PTMGPT2-parallel.py` on line 39 of the `Snakefile` (within the rule `predict_PTMs`) with `run-PTMGPT2.py`, which will look for a single GPU only. 
