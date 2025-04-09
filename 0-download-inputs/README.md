### Download input files used by downstream pipeline steps

This section of the pipeline handles the download of required files to `data-files/`:

* orf_trans.fasta - SGD yeast open reading frames
* UP000002311_559292_YEAST_v4.tar - AlphaFold2 structure predictions for the yeast proteome
* esm_if1_gvp4_t16_142M_UR50.pt 0 ESM-IF model used by Cagiada predictor in 5-dG-calculations
* YEAST_559292_idmapping.dat - Uniprot ID mapping file for yeast
* uniprot_sprot.dat - UniProt database 2025_01

After downloading the files, the second rule will extract SEQRES records from each AF2 PDB file and convert them to fasta format

To run this pipeline step in isolation, run the command `snakemake` from the `0-download-inputs` directory.

*N.B.*: This pipeline will add approximately 7.3 GB of data to the `data-files/` directory.
