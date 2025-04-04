### Parse UniProt database and add GO terms and additional annotation to each protein

Read in the file `uniprot_sprot.dat` and parse it to extract:

* Subcellular location data
* Function annotations
* Post-translational modifications

Execute this portion of the pipeline by running the command

`snakemake --use-conda`

Note: This code should be considered experimental; it requires additional manual checking by comparison to UniProt database entries

Snakemake pipeline currently fails during add-uniprot-info.py step, something weird going on to be debugged.
