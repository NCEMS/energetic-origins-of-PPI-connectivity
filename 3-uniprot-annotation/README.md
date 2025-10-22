### `3-uniprot-annotation`: Parse UniProt database and add GO terms and additional annotation to each protein

Corresponding configuration file section: `ANNOTATION WITH UNIPROT INFO`

This pipeline will:

1. Extract records from `uniprot_sprot.xml` for your organism of choice
2. Process these records to get subcellular location, function annotation, post-translational modifications, and GO terms
3. Optionally, process and integrate ProteomeXchange post-translational modification data (if a file is not provided, empty columns will be inserted into the output pd.DataFrame)
4. Output an updated annotated node pd.DataFrame

The UniProt search and ProteomeXchange data used are controlled by `../config-files/s288c.config`

N.B., you must supply the exact string used within UniProt to identify your organism of interest. For S. cerevisiae, this is "Saccharomyces cerevisiae (strain ATCC 204508 / S288c)". If you get an empty `{ORGANISM_TAG}-uniprot_sprot-parsed.csv` file, confirm your organism string is correct and try again. 
