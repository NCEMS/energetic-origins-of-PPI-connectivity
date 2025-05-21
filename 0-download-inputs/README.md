### Download input files used by downstream pipeline steps

This section of the pipeline handles the download of required files to `data-files/`:

* orf_trans.fasta - SGD yeast open reading frames
* UP000002311_559292_YEAST_v4.tar - AlphaFold2 structure predictions for the yeast proteome
* esm_if1_gvp4_t16_142M_UR50.pt - ESM-IF model used by Cagiada predictor in 5-dG-calculations
* YEAST_559292_idmapping.dat - Uniprot ID mapping file for yeast
* uniprot_sprot.xml - UniProt database 2025_01
* go.obo - Gene Ontology term library
* Yeast_GSB_phospho_all_prots_0125.csv - PTMeXchange phosphorylation site information for yeast

After downloading the files, the second rule will extract SEQRES records from each AF2 PDB file and convert them to fasta format

This step will require 60-90 min depending on connection speeds and the write speed of your file system. 

*N.B.*: This step will add approximately 7.3 GB of data to the `data-files/` directory.


