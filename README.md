
![](images/center-sized.jpg)
![](images/NSF-NCEMS-Blue.svg)

## Background and structure

### Repository contents

* This repository contains the code required to create, analyze, and explore the ANotated Yeast Interactome (ANYI), a heavily annotated yeast protein-protein interaction dataset. 
* Steps zero through eighteen (each with a correspond directory from `0-download-inputs` through to `18-finalize`) are run in series to produce the annotated interactome, saved in the file `18-finalize/processed-data/20260210-s288c-ANnnotated-Yeast-Interactome.pkl`
* In addition to the code required to generate this file, we also include a Docker image and interactive browser tool, ANYI Browser. See the Quick Start instructions below for details.

### Repository structure

* Nineteen directories with a numerical prefix {0, ..., 18} constitute the 19 steps required to assembled the ANotated Yeast Interactome. Each contains an individual Snakemake pipeline that carries out an action like downloading input data or adding an annotation to the base node file. 
* The directory `figures` contains sub directories with Jupyter Notebooks that generate the figures for [INSERT CITATION]
* The directory `docker` contains code and examples associated with the ANYI Browser tool

### Ways to use this reposistory

* There are three ways of using this repository and its associated Docker image:
	* Run Mode 1 (Quick start) - you can clone this repository and the Docker image onto your machine following the `Quick start` instructions below to launch the ANYI Browser tool to visualize yeast PPIs. 
	* Run Mode 2 (Reproduce key results) - you can use the repository and Docker image to reproduce all figures and key results without dealing with licensing agreements and expensive calculations
	* Run Mode 3 (Complete pipeline run) - you can rerun the entire pipeline including expensive calculations

### Files and directories

* The repository root directory contains the following files and folders:

**Table 1. Repository directories and files**
| Step Number | Name | Description |
|-----:|:----:|:------------|
|1|0-download-inputs| Download, unpack, and pre-process inputs |
|2|1-network-centrality| Compute network centrality metrics |
|3|2-sequence-parsing| Add sequence information, predict transmembrane proteins, predict signal sequences |
|4|3-uniprot-annotation| Add UniProt localization, post-translational modification, etc. data |
|5|4-idr-properties| Predict disordered regions and their sequence and dynamical properties |
|6|5-dG-calculations| Predict dG for each protein with empirical model and ESM-IF model |
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
|N/A|config-files| Contains the configuration file used to control inputs and outputs for Snakemake |
|N/A|docker| Contains information needed to build the Docker container associated with this repository |
|N/A|figures| Contains subdirectories corresponding to all files in [MANUSCRIPT LINK] |
|N/A|images| Contains a few images rendered in the README.md |
|N/A|`.dockerignore`| Files and folders not to be included in the Docker build |
|N/A|`.gitignore`| Files and folders not to be included in the Git repo |
|N/A|`README.md`| The file you are reading now |
|N/A|`docker_build.sh`| Contains the bash command used to build the Docker image |
|N/A|`full-rerun.sh`| Bash script that automates rerunning the entire pipeline |
|N/A|`minimal-rerun.sh`| Bash script that automates rerunning the entire pipeline except for expensive steps / steps requiring licensed software |

* All folders include their own `README.md` files explaining their contents. 

### Note on runtimes and requirments

* All steps of this pipeline were executed on an Ubuntu 22.04 machine with:
	* 2 x NVIDIA RTX 6000 Ada Gene GPUs
	* 112 threads on 56 Intel(R) Xeon(R) w9-3495X CPUs
	* 1 TB of memory
* Runtimes in the `README.md` files of individual pipeline steps, e.g. `1-netowork-centrality/README.md`, refer to the expected runtime on an equivalent system. 
* All input data requires ~110 GB of storage; all intermediate files and outputs bring the total size to ~250 GB.

## Run Mode 1 (Quick Start)

* This run option lets you immediately start interacting with ANYI

### Using the ANotated Yeast Interactome (ANYI) Browser tool

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
docker pull dannissleypsu/anyi-browser:v1.0.0
```

**Step 3** - Launch JupyterLab

* Run the command below to launch JupyterLab in the environment required by the ANYI Browser tool.

```bash
docker run --rm -it -p 8888:8888 \
  -e NB_UID=$(id -u) -e NB_GID=$(id -g) \
  -v "$PWD":/home/jovyan/work \
  dannissleypsu/anyi-browser:v1.0.0
```

* Once you have run the command above, copy the URL from your terminal into a web browser window.
* You can then use the navigation pane on the left to enter the `docker` folder and then `ANYI-browser` and then open `ANYI-browser.ipynb`. 
* By executing the code cells in this notebook and then clicking the `Launch` button, you can interact with the annotations in ANYI as well as their protein structures and key proteostasis metrics.

## Run Mode 2 (Reproduce Key Results)

* This run option allows you to skip dealing with licensed software and long runtimes by using precomputed data while still following the assembly process for ANYI.
* Following these instructions will allow for the generation of the same data currently in `18-finalize/processed-data/20260210-s288c-ANnnotated-Yeast-Interactome.pkl`
* In some places, you will need to update `config-files/minimal-rerun-s288c.config` while following the below instructions.

**Step 1 - Clone this repository**
```bash
git clone https://github.com/NCEMS/energetic-origins-of-PPI-connectivity.git
cd energetic-origins-of-PPI-connectivity
```

**Step 2 - Create the environment**

* You can also use the default environment from Run Mode 1 (i.e., within the Docker container) to execute this run mode. If you would prefer to build your own environment, however, use the below commands. 

```bash
conda env create -f docker/environment.yaml
conda activate anyi
```

**Step 3 - Setup `gocommands`**

* Follow the instructions to install `gocommands` in `0-download-inputs/README.md`
* Insert the absolute path to the `gocmd` executable into the file `config-files/minimal-rerun-s288c.config` for the term `gocommands` so that you have the line `gocommands: "/absolute/path/gocmd"` in the `download_inputs` section
* This series of steps enables the Snakemake pipeline in `0-download-inputs` to pull data from the CyVerse data store

**Step 4 - Download required data from CyVerse**

* To avoid rerunning expensive calculations we need to download some additional data from CyVerse. These downloads represent outputs from particularly long calculations.

(1) Rosetta data: from the repo root directory, run the following commands:

```bash
cd 6-Rosetta-scoring
gocmd get --progress /iplant/home/shared/NCEMS/working-groups/energetic-origins/required-data/6-Rosetta-scoring/Rosetta-N10_20251016.tar.gz .
tar -xvf Rosetta-N10_20251016.tar.gz
touch processed-data/scores/.all_scores_done
```

(2) PTMGPT2 data: from the repo root directory, run the following commands:

```bash
cd 10-predict-PTMs
gocmd get --progress /iplant/home/shared/NCEMS/working-groups/energetic-origins/required-data/11-predict-PTMs/PTMGPT2-predictions.tar.gz .
tar -xvf PTMGPT2-predictions.tar.gz
```

(3) Updated AlphaFold2 predictions needed by the pipeline: from any directory not in the repository, run the commands:

```bash
gocmd get --progress /iplant/home/shared/NCEMS/working-groups/energetic-origins/required-data/5-dG-calculations/alphafold2-structures .
```

After this download completes, update the `AF2_dir` parameter in `config-files/rerun-s288c.config` to point to the `alphafold2-structures/s288c` directory in your downloads folder, like this:

```json
AF2_dir: "/absolute/path/alphafold2_structures/s288c"
```

**Step 5 - Run the pipeline**

```bash
bash minimal-rerun.sh config-files/minimal-rerun-s288c.config
```

* This helper script will rerun all pipeline steps except for slow calculations or calculations requiring licensed software in `2-sequence-parsing`, `5-dG-calculations`, `6-Rosetta-scoring`, and `10-predict-PTMs`.
* The final output will be written to `18-finalize/processed-data/minimal-rerun-20260210-s288c-ANnnotated-Yeast-Interactome.pkl`
* The `20260210-s288c` part of this file name indicates the current ANYI build for the database
* This run will take ~60 min, most of which is required to download data

## Run Mode 3 (Complete Pipeline Run)

* This run option will require significant time and computational resources (approximately 45 days with 112 CPUs)
* Following these instructions will generate slightly different results from the data currently in `18-finalize/processed-data/20260210-s288c-ANnnotated-Yeast-Interactome.pkl` for Rosetta calculations
* In some places, you will need to update `config-files/full-rerun-s288c.config` while following the instructions below. 

**Important note**: the Snakemake pipeline in `0-download-inputs/Snakefile` will automatically download fixed versions of certain input data like, for example, open reading frame sequences from SGD. If you want to use updated versions, you will need to manually download them and then update `config-files/full-rerun-s288c.config` to point to these new files. 

**Step 1 - Clone this repository**

```bash
git clone https://github.com/NCEMS/energetic-origins-of-PPI-connectivity.git
cd energetic-origins-of-PPI-connectivity
```

**Step 2 - Create the environment**

```bash
conda env create -f docker/environment.yaml
conda activate anyi
```

**Step 3 - Setup `gocommands`**

* Follow the instructions to install `gocommands` in `0-download-inputs/README.md`
* Insert the absolute path to the `gocmd` executable into the file `config-files/rerun-s288c.config` for the term `gocommands` so that you have the line `gocommands: "/absolute/path/gocmd"` in the `download_inputs` section

**Step 4 - Download required software**

1) Rosetta

* Academic users can download Rosetta for free [here](https://downloads.rosettacommons.org/software/academic/)
* Once you have Rosetta installed, insert the absolute path to `relax.static.linuxgccrelease` into `config-files/full-rerun-s288c.config` in the `Rosetta_scoring` section:

```json
relax_executable: "/path/to/rosetta.binary.ubuntu.release-371/main/source/bin/relax.static.linuxgccrelease"
```

2) SignalP 6.0

* Academic users can navigate to the webpage [here](https://services.healthtech.dtu.dk/services/SignalP-6.0/) and sign the license agreement to be sent download information. SignalP Version 6.0h fast was used for the calculations included in this repository.
* Once you have downloaded SignalP and unpacked the download, you will have a directory named `signalp6_fast` or similar. Move this directory in `2-sequence-parsing/python-scripts/` so that, from the repo root directory, you can `ls 2-sequence-parsing/python-scripts/signalp6_fast`
* Now, copy the weights into the expected directory:

```bash
cp 2-sequence-parsing/python-scripts/signalp6_fast/signalp-6-package/models/distilled_model_signalp6.pt 2-sequence-parsing/python-scripts/signalp6_fast/signalp-6-package/signalp/model_weights/
```

Finally, edit `2-sequence-parsing/env/signalp-fast-local.yml` such that this section:

```json
  - pip:
      - -e /absolute/path/to/energetic-origins-of-PPI-connectivity/2-sequence-parsing/python-scripts/signalp6_fast/signalp-6-package
```

includes an absolute path to the `signalp-6-package` folder you just copied into your `2-sequence-parsing/python-scripts/` directory. 

3) PTMGPT2

* To enable PTMGPT2 runs, we need to clone the repository and download the models from Zenodo and then direct our configuration file to them. 
* In a directory outside of this repository, run the command:

```bash
git clone 
https://github.com/pallucs/PTMGPT2.git
```

This repository includes the `Tokenizer` folder. Next, download the required models (again, do this outside of this repository):

```bash
mkdir PTMGPT-2-models
cd PTMGPT-2-models
wget https://zenodo.org/records/11371883/files/PTMGPT2-models-Part1.zip
wget https://zenodo.org/records/11362322/files/PTMGPT2-models-Part2.zip
gunzip https://zenodo.org/records/11371883/files/PTMGPT2-models-Part1.zip
gunzip https://zenodo.org/records/11362322/files/PTMGPT2-models-Part2.zip
```

Then, update `config-files/full-rerun-s288c.config` to have the lines:

```json
gpt_model_path: "/absolute/path/to/PTMGPT-2-models"
tokenizer_path: "/absolute/path/to/PTMGPT2/Tokenizer"
```

Note that we direct `gpt_model_path` to a single directory containing all of the contents of `PTMGPT2-models-Part1.zip` and `PTMGPT2-models-Part2.zip` named `PTMGPT-2-models`

4) AlphaFold2 structure predictions

* By default, the pipeline compares reference sequences from the SGD to the sequences in the supplied AlphaFold2 predictions and, if there is a mismatch, tries to find a "rescue" structure in the path provided at the `AF2_dir` variable in the configuration file. If you want to rerun this step as well, you will need to supply AlphaFold2 structures matching the proteins in:

`/iplant/home/shared/NCEMS/working-groups/energetic-origins/required-data/5-dG-calculations/alphafold2-structures`
* You will then need to supply the absolute path to the directory containing protein-level directories each containing predictions in `config-files/full-rerun-s288c.config` for the variable `AF2_dir`.
* If you want to omit only this step, follow the instructions in `(3) Updated AlphaFold2 predictions needed by the pipeline` at the bottom of the `Run Mode 2` section of the README above.

**Step 5 - execute the complute pipeline**

* You can run all required commands in series by executing the helper bash script `full-rerun.sh`:

```bash
bash full-rerun.sh
```

## License and Attribution

This repository is an output of the National Synthesis Center for Emergence in the Molecular and Cellular Sciences (NCEMS) and is associated with the **Energetic Origins of Connectivity Within Protein Interaction Networks Working Group**.

Unless otherwise noted, the source code in this repository is licensed under the **GNU General Public License, version 3 or later (GPL-3.0-or-later)**. Please refer to the `LICENSE` file for the complete license text.

Users of this repository should provide appropriate attribution to **NCEMS** and the **Energetic Origins of Connectivity Within Protein Interaction Networks Working Group** in derivative works, presentations, publications, and other reuse, where applicable.

Please appropriately cite the relevant manuscript(s) arising from this Working Group:

- Nissley DA, Goel M, Castellanos-Girouard X, Kuntz CP, Wang Y, Mukhtar MS, Serohijos A, Schlebach JP. **ANYI: The ANnotated Yeast Interactome**. *Manuscript in preparation*.
- Goel M, Nissley DA, Castellanos-Girouard X, Kuntz CP, Wang Y, Mukhtar MS, Serohijos A, Schlebach JP. **Protein Stability, Turnover Kinetics, and Abundance Constrain the Scaling of Protein Interaction Networks**. *Manuscript in preparation*.
