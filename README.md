
### SETUP

1. You will need to first update `conda` on CyVerse by running the command:

`conda update -n base -c defaults conda`

(CyVerse appears to have 24.3.0, we need 24.7.1 or later for Snakemake to work correctly)

2. Setup a conda environment with Snakemake by running the command

`conda create --name snakemake -c bioconda -y snakemake`

followed by the command

`conda activate snakemake`

3. You can now run the pipeline by entering the command `./run_pipeline.sh <.config file>`

The `.config` file contains all commonly changed parameters, including those used to label output files. 

### PIPELINE

You can run the pipeline with the command:

`./run_pipeline.sh [.config file]`

For example, 

`.run_pipeline.sh config-files/s288c.config`

#### Step 0 - Download inputs

Downloads the required input data (e.g., fasta protein sequences, ESM-IF weights, etc.) and extracts AlphaFold2 structure sequences. 

#### Step 1 - Calculate network centrality metrics for each node

Uses NetworkX to annotate network with centrality metrics

#### Step 2 - Add sequence information for each node

Adds sequence information as possible to each node as well as DeepTMHMM annotations (predicts if proteins are TM, secreted, globular, etc.)

#### Step 3 - Add UniProt & GO annotations

Parses the UniProt database and inserts annotation information on function, subcellular location, post-translational modifications, and gene ontology terms

#### Step 4 - Predict IDRs and their properties

Predicts IDRs (metapredict v3.0) and annotates each of them with sequence- (CIDER) and ensemble-based (ALBATROSS) parameters

#### Step 5 - Carry out dG calculations

Uses Eq. 1 from Ghosh & Dill 2010 to predict dG for each protein sequence (minus cleaved signal sequences)

If requested in the .config file, will also run Cagiada stability predictions (required ~4.5 h for entire Yeast Interactome dataset; uses GPU)

#### Step 6 - Flatten the database

Converts the database from "one row per protein" to "one row per IDR"
