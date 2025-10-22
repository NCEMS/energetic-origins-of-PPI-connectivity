### `2-sequence-parsing`: Add sequence information and DeepTMHMM annotations to nodes

Corresponding configuration file section: `SEQUENCE PARSING`

This pipeline will:

(1) Add sequence information for genes as possible (sourced from SGD; see `../0-download-inputs/README.md` for details)
(2) Add DeepTMHMM predictions for each sequence (DeepTMHMM currently must be run outside of the pipeline and the path to the results included in the .config)
(3) Predict signal peptide cleavage sites with SignalP
(4) Add updated sequence information in "trimmed_sequence" column of the annotated node network

SignalP, when able to access a GPU, will require ~15 min. All other rules in this step require <2 min. 

To run SignalP after downloading and unpacking the directory as described in `../README.md`, follow these steps:

* You *must* copy the model weights from the model folder into the signalP folder:

`cp 2-sequence-parsing/python-scripts/signalp6_fast/signalp-6-package/models/distilled_model_signalp6.pt 2-sequence-parsing/python-scripts/signalp6_fast/signalp-6-package/signalp/model_weights/`

* You *must* edit the environment file `env/signalp_fast_local.yml` to insert the absolute path to the directory `signalp-6-package`. In the version distributed with this repository, the path is `/home/dan182/energetic-origins/energetic-origins-of-PPI-connectivity/2-sequence-parsing/python-scripts/signalp6_fast/signalp-6-package`

DeepTMHMM notes:

* The directory `DeepTMHMM-runs` contains the pre-computed results for the S288C yeast open reading frames.
* The file `DeepTMHMM-runs/command` provides the command used to generate the results in `DeepTMHMM-runs/s288c-results/all-predictions-s288c.3line`
* To rerun this section, install DockerHub and the image `dtu/deeptmhmm:1.0.24`. Start the container and then run the command in `DeepTMHMM-runs/command`
* Processing all sequences will require 24-36 hours without a GPU. 
