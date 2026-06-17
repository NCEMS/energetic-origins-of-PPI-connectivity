### `10-translation-efficiency`: Incorporate translation efficiency data

Corresponding configufle section: `Translation efficiency`

Currently getting CondaHTTPError when trying to build the environment:

```text
Assuming unrestricted shared filesystem usage.
host: E1-055815-ncems
Building DAG of jobs...
Creating conda environment 10-translation-efficiency/env/add-translation-speed.yml...
Downloading and installing remote packages.
CreateCondaEnvironmentException:
Could not create conda environment from /home/dan182/energetic-origins/arabidopsis-branch/energetic-origins-of-PPI-connectivity/10-translation-efficiency/env/add-translation-speed.yml:
Command:
conda env create --quiet --no-default-packages --file "/home/dan182/energetic-origins/arabidopsis-branch/energetic-origins-of-PPI-connectivity/.snakemake/conda/a2622a8aa26b4455c7b7408bcf3354cf_.yaml" --prefix "/home/dan182/energetic-origins/arabidopsis-branch/energetic-origins-of-PPI-connectivity/.snakemake/conda/a2622a8aa26b4455c7b7408bcf3354cf_"
Output:
Collecting package metadata (repodata.json): ...working... failed

CondaHTTPError: HTTP 429 TOO MANY REQUESTS for url <https://conda.anaconda.org/conda-forge/linux-64/repodata.json>
Elapsed: 00:06.197027
CF-RAY: a0d4d50e0e2572ed-IAD

An HTTP error occurred when trying to retrieve this URL.
HTTP errors are often intermittent, and a simple retry will get you on your way.
'https//conda.anaconda.org/conda-forge/linux-64'
```

### Legacy yeast information

Anticipated execution time: 1 min

To run this pipeline step in isolation, run the command:

```bash
snakemake -c all --use-conda --conda-frontend conda --snakefile 9-translation-efficiency/Snakefile --configfile config-files/s288c.config
```


* Merges translation efficiency information from `scikit-ribo` into the dataset
* Original data are contained within the `scikit-ribo_manuscript` repo [here](https://github.com/schatzlab/scikit-ribo_manuscript/blob/master/Data/skr_weinberg_genesTE.csv). 
