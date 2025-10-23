### `18-Y2H-data`: Add Y2H-derived PPI network to the existing network

Corresponding configuration file section: `Get Y2H network info & integrate`

The source file from CSBB `../0-download-inputs/data-files/Y2H_union.txt` was downloaded from `https://interactome.dfci.harvard.edu/S_cerevisiae/download/Y2H_union.txt`

This file was cleaned using the below steps:

1. Create the required environment: `conda env create -f env/clean-Y2H.yml`
2. Activate the environment: `conda activate clean-Y2H`
3. Run the command `python python-scripts/clean-Y2H-union.py --input_edges ../0-download-inputs/data-files/Y2H_union.txt --output_edges ../0-download-inputs/data-files/Y2H_union-clean.txt`

The file `../0-download-inputs/data-files/Y2H_union-clean.txt` was then loaded into Cytoscape to (1) run the `_wkshell` app, (2) create a node list, and (3) create an edge list

The resulting files are `../0-download-inputs/data-files/Y2H_union-clean_nodes.txt` and `../0-download-inputs/data-files/Y2H_union-clean_edges.txt`; once these files are in hand, you are ready to run this pipeline step (none of this is necessary in practice, as the node and edge files are distributed within the repo). 

Note that this section of the pipeline uses scripts in `../1-network-centrality`; the paths in your .config file should be relative paths from this directory (`18-Y2H-data`) to the directory `../1-network-centrality/python-scripts/...`

Note that the Y2H network is merged into the Yeast Interactome network to enable comparisons between nodes appearing in both datasets. However, this merge step means that *the Y2H data incorporated into the Yeast Interactome are incomplete*, as nodes that appear in the Y2H network but not the Yeast Interactome will be omitted. 
