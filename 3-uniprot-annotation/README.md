### Parse UniProt database and add GO terms and additional annotation to each protein

This pipeline will:

1. Extract records from `uniprot_sprot.xml` for your organism of choice
2. Process these records to get subcellular location, function annotation, post-translational modifications, and GO terms
3. Optionally, process and integrate ProteomeXchange post-translational modification data
4. Output an updated annotated node pd.DataFrame

The UniProt search and ProteomeXchange data used are controlled by `../config-files/s288c.config`
