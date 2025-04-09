### Parse UniProt database and add GO terms and additional annotation to each protein

Read in the file `uniprot_sprot.dat` and parse it to extract:

* Subcellular location data
* Function annotations
* Post-translational modifications

To run this pipeline step in isolation, run the command `snakemake --use-conda` from the `3-uniprot-annotation` directory.

*N.B.*: This code should be considered experimental; it requires additional manual checking by comparison to UniProt database entries
