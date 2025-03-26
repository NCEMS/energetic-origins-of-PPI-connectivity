### DESCRIPTION

Code and some data files for annotating a network with information about protein stability and network centrality metrics. 

The file data-files/The_Yeast_Interactome.cys was downloaded from [here](http://yeast-interactome.biochem.mpg.de:3838/interactome/session/f06cae8f2dbb9777aacde99863e20b61/download/downloadCyto?w=) and converted to edge and node graphs (as found in data-files/) using Cytoscape to export the data in format readable by the Python module NetworkX. 

### SETUP

To run this code on CyVerse, follow the steps below. 

* Log into CyVerse and navigate to the [Discovery Environment](de.cyverse.org)
* On the Apps tab, select "Jupyter Lab PyTorch GPU". An instance with at least 8 CPUs and 64 GB of memory alongside the GPU is suggested.
* Once you instance has launched, open a Terminal and navigate to the directory `/home/jovyan/data-store/` and run the command `git clone https://github.com/NCEMS/energetic-origins-of-PPI-connectivity.git` to clone the repo. (Note well - you will need to setup an SSH key to enable pulling code from this private repo).
* Once the code is downloaded, navigate into the repo directory and run the command `./bash-scripts/env-setup.sh` to build the necessary environments for the pipeline
* After the environments have been constructured, run the command `conda activate snakemake` followed by the command `./snakemake.command` to execute the pipeline
  * You can examine the file `Snakefile` to view a summary of the steps performed during the run. 

### USAGE

This pipeline uses Snakemake to integrate pipeline steps, currently amounts to:

The current version of this pipeline is designed to run in the "Jupyter Lab PyTorch GPU" app on CyVerse. The "cagiada_stability" section will give an error if you attempt to run it without an available GPU using CUDA. It can be adapted to run on CPUs, but is ~10x slower. 

0. Setup required environments on CyVerse
1. Download the required inputs (reference proteins sequences, esm model for Cagiada stability, yeast ID mappings from UniProt, etc.)
2. Use DeepTMHMM to annotate sequences with signal peptide/TM segment information
2. Load network and perform IDP annotation with metapredict
3. Add network centrality metric information
4. Run Cagiada stability as requested for proteins and add to network nodes

To do:

-1. Add typing checks to the code with the typing library
0. Compare P20484,YKL021C sequences between AF2 and SGD ref, appears to be source of error in cagiada_stability
1. Run cagiada_stability predictions for all proteins in nodes_df
2. Add DisProt-based cutoff for determining which proteins are/are not disordered
2. Add CentralityCosDist to list of network calculations. Python implementation appears to be available here: https://github.com/nilesh-iiita/CentralityCosDist
3. Need to handle ambiguous node names of the form "gene1;gene2;gene3" that appear in The_Yeast_Interactome files (currently left in but not annotated)

FIXED/DONE issues

* Having some trouble with metapredict install from environment.yml file; not sure exactly what is going on (fixed, just syntax problems) (FIXED, pip syntax issue)
* Add test for cagaida_stability.py; make sure prediction for test protein matches results from author's Google Colab
* Incorporate DeepTMHMM predictions into the pipeline 
* Add test that makes sure sequence in nodes_df is the same length as the sequence in the AF2 model used to predict stability

Known issues:

1. (Cagiada stability) Predictions on multi-domain proteins or proteins with complex folding kinetics show an absolute stability overestimated compared to the real one.
2. (Cagiada stability) Predictions are limited to proteins with 1023 residues (Max protein size for the ESM-IF language model)
