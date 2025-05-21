### Predict post-translation modifications

* Runs the model PTMGPT2 to predict the locations of 19 different post-translational modifications

* Predictions are run on SignalP-trimmed sequences, but residue indices are mapped back to the original numbering for the untruncated sequence. 

*Note*: The current Snakefile assumes you have 2 GPUs with device IDs [0, 1]. Using two RTX6000 GPUs in coarse-grain parallel, predictions required ~48 hours for all 19 models. 
