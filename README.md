Code and some data-files for network analyses and protein stability predictions. 

The file data-files/The_Yeast_Interactome.cys was downloaded from [this site](http://yeast-interactome.biochem.mpg.de:3838/interactome/session/f06cae8f2dbb9777aacde99863e20b61/download/downloadCyto?w=) and converted to edge and node graphs (as found in data-files/) using Cytoscape to export the data. 

This pipeline uses Snakemake to integrate pipeline steps, currently amounts to:

0. Setup required environments DONE (MacOS & CyVerse)
1. Download the required inputs (reference proteins sequences, esm model for Cagiada stability, yeast ID mappings from UniProt, etc.) DONE
2. Use DeepTMHMM to annotate sequences with signal peptide/TM segment information (TO-DO, next on the list)
2. Load network and perform IDP annotation with metapredict DONE
3. Add network centrality metric information DONE
4. Run Cagiada stability as requested for proteins and add to network DONE

To be added:

1. Incorporate DeepTMHMM calculations to identify transmembrane and secreted proteins. Transmembrane proteins should probably not have their disorder scores calculated, and signal peptides should be omitted from disorder calculations (per email from Jon S. on March 17, 2025)
	This requires a Docker container; will need to get this running using the nvidia GPUs on the "Jupyter_Lab_PyTorch_GPU_analysis1" App of CyVerse for efficiency
	Not clear if I have sudo on Cyverse to install docker, may need to build container on my local machine (annoying...)
2. Incorporate test that makes sure sequence in nodes_df is the same length as the sequence in the AF2 model used to predict stability (have FASTA files, just need to add the check)
3. Need to handle ambiguous node names of the form "gene1;gene2;gene3" that appear in The_Yeast_Interactome files

FIXED/DONE issues 

* Having some trouble with metapredict install from environment.yml file; not sure exactly what is going on (fixed, just syntax problems) (FIXED, pip syntax issue)

Know issues:

1. (Cagiada stability) Predictions on multi-domain proteins or proteins with complex folding kinetics show an absolute stability overestimated compared to the real one.
2. (Cagiada stability) Predictions are limited to proteins with 1023 residues (Max protein size for the ESM-IF language model)
