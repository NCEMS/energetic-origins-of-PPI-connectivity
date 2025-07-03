### Add Y2H-derived PPI network to the existing network

Note well: this piece of the pipeline breaks the automated flow; you must have already cleaned (see below) the input data and run it through Cytoscape to get `_wkshell` before running this step

The source file from CSBB `../0-download-inputs/data-files/Y2H_union.txt` was downloaded from https://interactome.dfci.harvard.edu/S_cerevisiae/download/Y2H_union.txt

This file was cleaned using the below steps:

1. Create the required environment: `conda env create -f env/clean-Y2H.yml`
2. Activate the environment: `conda activate clean-Y2H`
3. Run the command `python python-scripts/clean-Y2H-union.py --input_edges ../0-download-inputs/data-files/Y2H_union.txt --output_edges ../0-download-inputs/data-files/Y2H_union-clean.txt`

The file `../0-download-inputs/data-files/Y2H_union-clean.txt` was then loaded into Cytoscape to (1) run the `_wkshell` app, (2) create a node list, and (3) create an edge list

The resulting files are `../0-download-inputs/data-files/Y2H_union-clean_nodes.txt` and `../0-download-inputs/data-files/Y2H_union-clean_edges.txt`; once these files are in hand, you are ready to run this pipeline step (none of this is necessary in practice, as the node and edge files are distributed within the repo). 

Note that this section of the pipeline uses scripts in 1-network-centrality; the paths in your .config file should be relative paths from this directory (`18-Y2H-data`) to the directory `../1-network-centrality/python-scripts/...`

