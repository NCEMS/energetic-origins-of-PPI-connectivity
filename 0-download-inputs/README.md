### Download input files used by downstream pipeline steps

To setup a conda environment with Snakemake, run the command

`conda create --name snakemake -c bioconda -y snakemake`

followed by the command

`conda activate snakemake`

You can then run the pipeline by entering the command `snakemake`

This pipeline handles the downloading of required files to `data-files/`:

* orf_trans.fasta - SGD yeast open reading frames
* UP000002311_559292_YEAST_v4.tar - AlphaFold2 structure predictions for the yeast proteome
* esm_if1_gvp4_t16_142M_UR50.pt 0 ESM-IF model used by Cagiada predictor in 5-dG-calculations
* YEAST_559292_idmapping.dat - Uniprot ID mapping file for yeast
* uniprot_sprot.dat - UniProt database 2025_01

*Note*: This pipeline will add approximately 7.3 GB of data to the `data-files/` directory.
