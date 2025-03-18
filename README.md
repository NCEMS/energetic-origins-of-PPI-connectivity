Code and some data-files for network analyses and protein stability predictions. 

Uses Snakemake to integrate pipeline steps, currently amounts to:

1. Download the required inputs (reference proteins sequences, esm model for Cagiada stability, yeast ID mappings from UniProt, etc.)
2. Load network and perform IDP annotation with metapredict
3. Add network centrality metric information
4. Run Cagiada stability as requested for proteins and add to network

To be added:

1. Incorporate DeepTMHMM calculations to identify transmembrane and secreted proteins. Transmembrane proteins should probably not have their disorder scores calculated, and signal peptides should be omitted from disorder calculations (per email from Jon S. on March 17, 2025)
