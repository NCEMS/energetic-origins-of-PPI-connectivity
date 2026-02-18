### Background and structure

#### Repository structure

* Nineteen directories with a numerical prefix {0, ..., 18} constitute the 19 steps required to assembled the ANotated Yeast Interactome. Each contains an individual Snakemake pipeline that carries out an action like downloading input data or adding an annotation to the base node file. 
* There are three ways of using this repository and its associated Docker image. 
	* First, you can clone this repository and the Docker image onto your machine and follow the `Quick start` instructions below to launch the ANYI Browser tool to visualize yeast PPIs. 
	* Second, you can use the repository and Docker image to reproduce all figures and key results by rerunning pipeline steps with pre-computed data (e.g., Rosetta relaxed structures)
	* Third, you can rerun the entire pipeline including computationally expensive steps like Rosetta relaxation
* Each of the latter two uses of the repository are described in detail below. 

##### Files and directories

**Table 1. Repository directories and files**
|Step Number| Name | Description |
|-----:|:----:|:------------|
|1|0-download-inputs| Download, unpack, and pre-process inputs |
|2|1-network-centrality| Compute network centrality metrics |
|3|2-sequence-parsing| Add sequence information, predict transmembrane proteins, predict signal sequences |
|4|3-uniprot-annotation| Add UniProt localization, post-translational modification, etc. data |
|5|4-idr-properties| Predict disordered regions and their sequence and dynamical properties |
|6|5-dG-calculations| Predict dG for each protein with empirical model and ESM-IF generative model |
|7|6-Rosetta-scoring| Relax AlphaFold2 structures with Rosetta and score |
|8|7-protein-half-life| Add protein half-life data |
|9|8-protein-expression| Add protein expression data |
|10|9-translation-speed| Add protein translation efficiency data from `scikit-ribo`|
|11|10-predict-PTMs| Predict PTMs with `PTMGPT2`|
|12|11-entanglement| Add entanglement data |
|13|12-chaperones| Add chaperone data |
|14|13-oligomers| Add oligomerization state/complex membership data from Complex Portal |
|15|14-domain-annotations| Add domain annotations from InterPro|
|16|15-essentiality| Add SGD protein essentiality information |
|17|16-Y2H-data| Add yeast two-hybrid data from Yu et al. 2008 |
|18|17-meltome-atlas| Add Meltome Atlas thermal stability data |
|19|18-finalize| Post-process the annotated node network for easy analysis |
|N/A|`README.md`| The file you are reading now |
|N/A|config-files| Contains the configuration file used to control inputs and outputs for Snakemake |
|N/A|docker| Contains information needed to build the Docker container associated with this repository |
|N/A|`docker_build.sh`| Contains the bash command used to build the Docker image |
|N/A|figures| Contains subdirectories corresponding to all files in [MANUSCRIPT LINK] |
|N/A|`run-pipeline.sh`| Bash script that automates rerunning the entire pipeline |
|N/A|`reproduce-results.sh`| Bash script that automates reproducing results without rerunning long calculations |

* All folders include their own `README.md` files explaining their contents. 

#### Runtimes

* All steps of this pipeline were executed on an Ubuntu 22.04 machine with 2 x NVIDIA RTX 6000 Ada Gene GPUs and 112 threads on 56 Intel(R) Xeon(R) w9-3495X CPUs. Runtimes in the README.md files of individual pipeline steps, e.g. `1-netowork-centrality/README.md` refer to the expected runtime on an equivalent system. 

### Quick start

#### Using the ANotated Yeast Interactome (ANYI) Browser tool

* This repository is designed to be used with a pre-built Docker image that contains the full ANYI runtime environment (JupyterLab + required Python packages). 

**Step 1** - Clone this repository

* Run the command below to clone this repository and then enter its root directory.

```bash
git clone https://github.com/<your-org>/energetic-origins-of-PPI-connectivity.git
cd energetic-origins-of-PPI-connectivity
```

**Step 2** - Pull the Docker image

* Run the command below to pull the Docker image

```bash
docker pull dannissleypsu/anyi:0.1.0
```

**Step 3** - Launch JupyterLab

* Run the command below to launch JupyterLab in the environment required by the ANYI Browser tool.

```bash
docker run --rm -it -p 8888:8888 \
  -e NB_UID=$(id -u) -e NB_GID=$(id -g) \
  -v "$PWD":/home/jovyan/work \
  dannissleypsu/anyi:0.1.0
```

* Once you have run the command above, copy the URL from your terminal into a web browser window.
* You can then use the navigation pane on the left to enter the `docker` folder and then `ANYI-browser` and then open `ANYI-browser.ipynb`. 
* By executing the code cells in this notebook and then clicking the `Launch` button, you can interact with the annotations in ANYI as well as their protein structures and key proteostasis metrics.

### Reproducing key outputs and figures

Simple command to use 
### Computational requirements, dependencies, and benchmarks

#### Requirements

This repository uses a series of Snakemake pipelines to assemble an annotated protein-protein interaction network. Each individual pipeline consists of discrete Python and bash processing steps (rules within Snakemake). 

To run the complete pipeline, you will need:

1. An internet connection to download obligate input files
2. ~250 GB of storage space (for all inputs and outputs)
3. A CUDA-enabled GPU for Cagiada et al. 2025 ESM-IF-based dG predictions & PTMGPT2 post-translational modification prediciton
4. 60-120 CPUs to enable Rosetta relaxation of protein structures in a reasonable timeframe

#### Dependencies

Nearly all dependency issues will be handled by Snakemake automatically by building conda environments based on the files in each pipeline's `env` subdirectory. However, if you want to rerun all pipeline steps you will need to download and install additional software. Converesly, if you want to skip some steps and use pre-generated data to save time, you will need to download it from CyVerse. 

If you want to rerun everything, follow the download instructions in the table below to setup SignalP, Rosetta, FoldX, and PTMGPT2.

| Step Number | Description | Instructions |
|------------:|:-----------:|:------------:| 
| 2 | SignalP6.0 for prediction of protein signal sequences | Download [here](https://services.healthtech.dtu.dk/cgi-bin/sw_request?software=signalp&version=6.0&packageversion=6.0h&platform=fast) and unpack `signalp-6.0h.fast.tar.gz`  into `2-sequence-parsing/python-scripts`. You should have the path `2-sequence-parsing/python-scripts/signalp6_fast/signalp-6-package/` available from the repo root directory. See `2-sequence-parsing/README.md` for additional setup steps.|
| 6 | Rosetta for structure relaxation and scoring | Download from [Rosetta Commons](https://rosettacommons.org/software/download/) and insert the absolute path to `relax.static.linuxgccrelease` or equivalent into the .config file in the Rosetta scoring section for the variable `relax_executabele`. |
|10 | PTMGPT2 models for post-translational modification prediction | The models [Part 1](https://zenodo.org/records/11371883) and [Part 2](https://zenodo.org/records/11362322) can be downloaded from Zenodo. Both .zip files should be unpacked into one directory and the absolute path to this directory inserted into the "predict post-translational modifications" section of the .config file for the variable `gpt_model_path`; you must also download the Tokenizer from https://github.com/pallucs/PTMGPT2 and add its path to your config file |

If you want to use existing data for yeast, follow the instructions in the table below to download it from CyVerse.

| Step Number | Description | Instructions |
|------------:|:-----------:|:------------:|
| 5           | Additional AlphaFold2 protein stuctures predictions | Download with 'gocommands' from the path `/iplant/home/shared/NCEMS/working-groups/energetic-origins/required-data/5-dG-calculations/alphafold2-structures` and update the `AF2_struc` path in your configuration file |
| 6           | Pre-computed Rosetta relaxed structures and scores | Download with `gocommands` from the path `/iplant/home/shared/NCEMS/working-groups/energetic-origins/required-data/6-Rosetta-scoring/Rosetta-N10_20251016.tar.gz` and place the contents of the .tar.gz archive in `6-Rosetta-scoring/processed-data/scores` |
|10           | Pre-computed PTMGPT2 predictions for yeast proteins | Download with `gocommands` from the path `/iplant/home/shared/NCEMS/working-groups/energetic-origins/required-data/11-predict-PTMs/PTMGPT2-predictions.tar.gz` and place the contents of the .tar.gz archive in `11-predict-PTMs/processed-data`|

Check the step-specific subdirectories for any additional setup instructions for SignalP, Rosetta, FoldX, and PTMGPT2.

##### Using `gocommands` to get data from the CyVerse Data Store

Visit [this page](https://learning.cyverse.org/ds/gocommands/installation) and follow the installation instructions for your system. Once installed, you should have the executable `gocmd` in the folder where you ran the installation command. 
With `gocmd` available in your system, you can download data required from CyVerse like so:

`cd 6-Rosetta-scoring`
`gocmd get --progress /iplant/home/shared/NCEMS/working-groups/energetic-origins/required-data/6-Rosetta-scoring/Rosetta-N10_20251016.tar.gz .`
`tar -xvf Rosetta-N10_20251016.tar.gz`


#### Benchmarks

Most pipeline steps include simple procedures like loading, cleaning, and merging datasets together. Some, however, require more heavy-duty computation. 

The main computational bottlenecks are:

| Step Number | Description | Time |
|------------:|:-----------:|:----:|
|           0 | Download, unpacking, and pre-processing of input data | Requires up to 2 hours depending on connection speeds and write speed of drive |
|           5 | dG prediction from structure with Cagiada et al. 2025 method | ~4 h on A16; ~25 min on RTX 6000 Ada Gene |
|           6 | Rosetta structure relxation & scoring | ~24 days with 96 Intel(R) Xeon(R) w9-3495X CPUs with N = 10 replicates per protein |
|          11 | Prediction of post-translational modifications with PTMGPT2 | ~48 h on 2 x RTX 6000 Ada Gene GPUs in coarse-grain parallel |
|          16 | Extract domain annotations from ~100 GB file | Depending on on your system, 10 min - 2 h |

As you will read below (see the section **Running the pipeline now** below), you can skip these expensive calculations if you just want to rerun the pipeline as-is. If you are running for a new organism/new proteins, these calculations are a one-time cost. 

### SETUP

The pipeline for assembling the final data product is in fact an ensemble of pipelines. You will need to install `conda` and then setup a top-level environment with Snakemake. 

1. Update your version of conda:

`conda update -n base -c defaults conda`

2. Setup a conda environment with Snakemake by running the command

`conda create --name snakemake -c bioconda -y snakemake`

followed by the command

`conda activate snakemake`

Note: This gave me issues on a new Ubuntu machine; if you have any problems, try `conda create -n snakemake -c conda-forge -c bioconda "python>3.11" "snakemake>=9,<10" biopython
`.

3. Download SignalP, Rosetta, FoldX, & PTMGPT2 models (see **Dependencies** section above)

4. You can now run the pipeline by entering the command `./run-pipeline.sh <.config file>`

The `.config` file contains all commonly changed parameters, including those used to label output files. The current config file to run all steps is `config-files/s288c.config`.

### Running the pipeline now

Once you have downloaded all required code and/or data, you can run the pipeline with the command:

`./run-pipeline.sh [.config file]`

For example, 

`./run-pipeline.sh config-files/s288c.config`

If you would like to run a specific pipeline step in isolation, you can use a command like:

`snakemake --snakefile /snakefile/path/Snakefile --configfile /configfile/path/config.config --use-conda --conda-frontend conda -c all`

in which you must replace `/snakefile/path/Snakefile` and `/configfile/path/config.config` with correct relative paths. For example, from the main repo working directory we could run pipeline step 8 in isolation with the command:

`snakemake --snakefile 8-protein-half-life/Snakefile --configfile config-files/s288c.config --use-conda --conda-frontend conda -c all`
