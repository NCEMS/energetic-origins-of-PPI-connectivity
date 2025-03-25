Code and some data-files for network analyses and protein stability predictions. 

The file data-files/The_Yeast_Interactome.cys was downloaded from [this site](http://yeast-interactome.biochem.mpg.de:3838/interactome/session/f06cae8f2dbb9777aacde99863e20b61/download/downloadCyto?w=) and converted to edge and node graphs (as found in data-files/) using Cytoscape to export the data. 

This pipeline uses Snakemake to integrate pipeline steps, currently amounts to:

The current version of this pipeline is designed to run in the "Jupyter Lab PyTorch GPU" app on CyVerse. The "cagiada_stability" section will give an error if you attempt to run it without an available GPU using CUDA. It can be adapted to run on CPUs, but is ~10x slower. 

0. Setup required environments on CyVerse
1. Download the required inputs (reference proteins sequences, esm model for Cagiada stability, yeast ID mappings from UniProt, etc.)
2. Use DeepTMHMM to annotate sequences with signal peptide/TM segment information
2. Load network and perform IDP annotation with metapredict
3. Add network centrality metric information
4. Run Cagiada stability as requested for proteins and add to network nodes

To be added:

1. Incorporate test that makes sure sequence in nodes_df is the same length as the sequence in the AF2 model used to predict stability (have FASTA files, just need to add the check)
2. Run cagiada_stability predictions for all proteins in nodes_df
3. Need to handle ambiguous node names of the form "gene1;gene2;gene3" that appear in The_Yeast_Interactome files

FIXED/DONE issues

* Having some trouble with metapredict install from environment.yml file; not sure exactly what is going on (fixed, just syntax problems) (FIXED, pip syntax issue)
* Add test for cagaida_stability.py; make sure prediction for test protein matches results from author's Google Colab

Know issues:

1. (Cagiada stability) Predictions on multi-domain proteins or proteins with complex folding kinetics show an absolute stability overestimated compared to the real one.
2. (Cagiada stability) Predictions are limited to proteins with 1023 residues (Max protein size for the ESM-IF language model)
