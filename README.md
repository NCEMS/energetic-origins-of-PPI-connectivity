
### SETUP

Running the complete pipeline requires significant CPU and GPU resources. While some steps can be run efficiently on CyVerse, this is not recommended for the complete pipeline. 

1. Update your version of conda:

`conda update -n base -c defaults conda`

(We need 24.7.1 or later for Snakemake to work correctly)

2. Setup a conda environment with Snakemake by running the command

`conda create --name snakemake -c bioconda -y snakemake`

followed by the command

`conda activate snakemake`

3. Download SignalP, Rosetta, FoldX, & PTMGPT2 models

SignalP6.0 fast can be downloaded from [this site](https://services.healthtech.dtu.dk/cgi-bin/sw_request?software=signalp&version=6.0&packageversion=6.0h&platform=fast) after accepting the academic licensing agreement. The contents of the downloaded `signalp-6.0h.fast.tar.gz` should be unpacked into `2-sequence-parsing/python-scripts` to enable the environment associated with the SignalP Snakemake rule to build correctly. For example, you should have the path `2-sequence-parsing/python-scripts/signalp6_fast/signalp-6-package/` available from the repo root directory. 

Rosetta can be downloaded from [Rosetta Commons](https://rosettacommons.org/software/download/) free of charge. Insert the absolute path to `relax.static.linuxgccrelease` or equivalent into the .config file in the Rosetta scoring section for the variable `relax_executable`.

FoldX can be [downloaded](https://foldxsuite.crg.eu/) after making an account and accepting the academic license agreement. Insert the absolute path to the pre-compiled binary in the FoldX section of the .config file for the variable `executable`.

PTMGPT2 models [Part 1](https://zenodo.org/records/11371883) and [Part 2](https://zenodo.org/records/11362322) can be downloaded from Zenodo. Both of these .zip files should be unpacked into one directory and the absolute path to this directory inserted into the "predict post-translational modifications" section of the .config file for the variable `gpt_model_path`

4. You can now run the pipeline by entering the command `./run_pipeline.sh <.config file>`

The `.config` file contains all commonly changed parameters, including those used to label output files. The current config file to run all steps is `config-files/s288c.config`.

Approximate timings for individual pipeline steps are listed in the README.md files within subdirectories for each step

5. Optional - if you want to rerun the pipeline from step 0, download additional data for `6-Rosetta-scoring` from:

`/iplant/home/shared/NCEMS/working-groups/energetic-origins/additional-data/6-Rosetta-scoring/scores`

As well as additional data for `7-FoldX-scoring` from:

`/iplant/home/shared/NCEMS/working-groups/energetic-origins/additional-data/7-FoldX-scoring/scores`

And place these directories in `6-Rosetta-scoring/processed-data/scores` and `7-FoldX-scoring/processed-data/scores`

Finally, place the data from:

`/iplant/home/shared/NCEMS/working-groups/energetic-origins/additional-data/11-predict-PTMs/` in `11-predict-PTMs/processed-data` and you are ready to go without needing to rerun expensive calculations.


### PIPELINE

You can run the pipeline with the command:

`./run_pipeline.sh [.config file]`

For example, 

`.run_pipeline.sh config-files/s288c.config`

#### 0-download-inputs

Downloads the required input data (e.g., protein ORF sequences, ESM-IF weights, etc.) and extracts AlphaFold2 structure sequences. 

Note that some required files are not downloaded at runtime but are distributed with this GitHub repo.

#### 1-network-centrality

Uses NetworkX to annotate network with centrality metrics; note that weighted k-shell calculation results are included in 0-download-inputs/data-files/The_Yeast_Interactom_nodes.csv.

#### 2-sequence-parsing

Adds sequence information as possible to each node 

Adds DeepTMHMM annotations (predicts if proteins are TM, secreted, globular, etc.; this step is precomputed using the DeepTMHMM Docker container)

Adds SignalP6.0 identification of cleavage sites for signal peptides

#### 3-uniprot-annotation

Parses the UniProt .xml database and inserts annotation information on function, subcellular location, post-translational modifications, and gene ontology terms.

PTMeXchange phosphorlyation sites are also integrated. 

#### 4-idr-properties

Predicts IDRs (metapredict v3.0) and annotates each of them with sequence- (CIDER) and ensemble-based (ALBATROSS) parameters using SPARROW. 

#### 5-dG-calculations

Uses Eq. 1 from Ghosh & Dill 2010 to predict dG for each protein sequence (minus cleaved signal sequences)

If requested in the .config file, will also run Cagiada stability predictions (required ~1 h for entire Yeast Interactome dataset on an RTX4500 GPU)

#### 6-Rosetta-scoring

Takes AF2 structures (with signal sequences cleaved as necessary) and runs Rosetta FastRelax on each of them, generating a single relaxed pose and Rosetta scoring function data

#### 7-FoldX-scoring

Uses FoldX scoring function to compute various energy parameters including a total stability for each of the poses generated by Rosetta FastRelax. 

#### 8-protein-half-life

Integrates protein half-life data from 10.1016/j.celrep.2014.10.065

#### 9-protein-expression

Integrates protein expression data from 10.1016/j.cels.2017.12.004

#### 10-translation-speed

Integrates translation efficiency information computed using scikit-ribo on Weinberg 2016 dataset

#### 11-predict-PTMs

Runs the model PTMGPT2 to predict post-translational modifications for each protein. Current list of predicted PTMs can be found by checking main() in 11-predict-PTMs/python-scripts/run-PTMGPT2-parallel.py

#### flatten

This is the final step in the pipeline; performs some minor cleanup of redundant columns and saves the database in two forms with two file types for each form:

.csv & .pkl files with "*-nodes-final-per-node.csv" ending contain the output annotated dataset on a per-node or per-protein basis (i.e., one row per node or protein)

.csv & .pkl files with "*-nodes-final-per-IDR.csv" ending contain the output annotated dataset flattend over IDRs, giving a file with one row per IDR rather than one row per node

These files can be found in `flatten/processed-data/`
