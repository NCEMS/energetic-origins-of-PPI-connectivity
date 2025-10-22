### `1-network-centrality`: Annotate network with centrality metrics

Corresponding configuration file section: `NETWORK CENTRALITY CALCULATIONS`

This pipeline step computes the following centrality metrics:

* Degree centrality
* Betweenness centrality
* Eigenvector centrality
* Closeness centrality
* Load centrality
* Page rank
* Information centrality (run on the largest connected subgraph)
* CentralityCosDist (using each of the above as part of the vector; default is to use all nodes as seed nodes)

Cytoscape version 3.10.3 was used to convert `The_Yeast_Interactome.cys` to edge and node tables. These two files are the key inputs to this step (see the configuration file)

Weighted k-shell decomposition was pre-computed in Cytoscape using The Yeast Interactome as the input; this information is stored in `../0-download-inputs/data-files/The_Yeast_Interactome_nodes.csv`

The app `wk-shell-decomposition` was downloaded from the [Cytoscape App Store](https://apps.cytoscape.org/apps/wkshelldecomposition). 

Two files will be written by the Snakemake pipeline in `1-network-centrality/processed-data`: `{output_prefix}-CentralityCosDist-input.csv` and `{output_prefix}-{organism_label}-nodes-centrality.csv`. 

The values of {output_prefix} and {organism_label} are taken from the configuration file at runtime. 

`{output_prefix}-CentralityCosDist-input.csv` is an intermediate file written in the specific format required for Centrality Cosine Distance calculations. 

`{output_prefix}-{organism_label}-nodes-centrality.csv` is the output file from this step containing the annotated node network. Each row corresponds to a node, with columns for centrality metrics added for each node. 

Centrality Cosine Distance calculations are carried out using code from the [Mukhtar Lab](https://github.com/nilesh-iiita/CentralityCosDist). This repository is distributed in `1-network-centrality/python-scripts/CentralityCosDist`. The directory `1-network-centrality/test` contains input data used to perform a test at runtime that the CentralityCosDist code is functioning as expected. 
