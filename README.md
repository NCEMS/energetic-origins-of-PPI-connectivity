Setup a conda environment with Snakemake by running the command

`conda create --name snakemake -c bioconda -y snakemake`

followed by the command

`conda activate snakemake`

You can then run the pipeline by entering the command `./run_pipeline.sh <.config file>`

The `.config` file contains all commonly changed parameters, including those used to label output files. 
