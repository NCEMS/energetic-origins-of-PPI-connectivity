### Parse UniProt database and add GO terms and additional annotation to each protein

Read in the file `uniprot_sprot.xml` and parse it to extract:

* Subcellular location data
* Function annotations
* Post-translational modifications
* GO terms

Read in the file `Yeast_GSB_phospho_all_prots_0125.csv` and parse it to extract:

* Yeast phosphorylation sites from mass spectrometry data on PTMeXchange
