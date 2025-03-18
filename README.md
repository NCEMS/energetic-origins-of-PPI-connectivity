Code and some data-files for network analyses and protein stability predictions. 

The file data-files/The_Yeast_Interactome.cys was downloaded from [this site](http://yeast-interactome.biochem.mpg.de:3838/interactome/session/f06cae8f2dbb9777aacde99863e20b61/download/downloadCyto?w=) and converted to edge and node graphs (as found in data-files/) using Cytoscape to export the data. 

This pipeline uses Snakemake to integrate pipeline steps, currently amounts to:

0. Setup required environments
1. Download the required inputs (reference proteins sequences, esm model for Cagiada stability, yeast ID mappings from UniProt, etc.)
2. Load network and perform IDP annotation with metapredict
3. Add network centrality metric information
4. Run Cagiada stability as requested for proteins and add to network

To be added:

1. Incorporate DeepTMHMM calculations to identify transmembrane and secreted proteins. Transmembrane proteins should probably not have their disorder scores calculated, and signal peptides should be omitted from disorder calculations (per email from Jon S. on March 17, 2025)
