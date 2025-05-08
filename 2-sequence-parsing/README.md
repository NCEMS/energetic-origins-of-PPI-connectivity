### Add sequence information and DeepTMHMM annotations to nodes

(1) Add sequence information for genes as possible
(2) Add DeepTMHMM predictions for each sequence (DeepTMHMM must be run outside of the pipeline and the path to the results included in the .config)
(3) Predict signal peptide cleavage sites with SignalP
(4) Add updated sequence information in "trimmed_sequence" column of the nodes_df
