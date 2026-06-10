### `1-network-centrality`: Annotate network with centrality metrics

Corresponding configuration file section: `NETWORK CENTRALITY CALCULATIONS`

Anticipated execution time: ~1 min

This pipeline step computes the following centrality metrics:

* Degree centrality
* Betweenness centrality
* Eigenvector centrality
* Closeness centrality
* Load centrality
* Page rank
* Information centrality (run on the largest connected subgraph)

Current A. thaliana centrality calculations:

| Input Node file | Input Edge file | Config file | Centrality Metric Output File |
|:---|:---|:---|:---|
| `../0-download-inputs/data-files/merged-nodes.csv` | `../0-download-inputs/data-files/merged-edges.csv` | `../config-files/union-athaliana.config` | `processed-data/union-athaliana-step1.csv` |
| `../0-download-inputs/data-files/in-vivo-nodes.csv` | `../0-download-inputs/data-files/in-vivo-edges.csv` | `../config-files/in-vivo-athaliana.config` | `processed-data/in-vivo-athaliana-step1.csv` |
| `../0-download-inputs/data-files/in-vitro-nodes.csv` | `../0-download-inputs/data-files/in-vitro-edges.csv` | `../config-files/in-vitro-athaliana.config` | `processed-data/in-vitro-athaliana-step1.csv` |
| `../0-download-inputs/data-files/intersection-nodes.csv` | `../0-download-inputs/data-files/intersection-edges-in_vivo.csv` | `../config-files/intersection-in-vivo-athaliana.config` | `processed-data/intersection-in-vivo-athaliana-step1.csv` |
| `../0-download-inputs/data-files/intersection-nodes.csv` | `../0-download-inputs/data-files/intersection-edges-in_vitro.csv` | `../config-files/intersection-in-vitro-athaliana.config` | `processed-data/intersection-in-vitro-athaliana-step1.csv` |


### Legacy yeast README.md information

Cytoscape version 3.10.3 was used to convert `The_Yeast_Interactome.cys` to edge and node tables. These two files are the key inputs to this step (see the configuration file)

Weighted k-shell decomposition was pre-computed in Cytoscape using The Yeast Interactome as the input; this information is stored in `../0-download-inputs/data-files/The_Yeast_Interactome_nodes.csv`

The app `wk-shell-decomposition` was downloaded from the [Cytoscape App Store](https://apps.cytoscape.org/apps/wkshelldecomposition). 

Two files will be written by the Snakemake pipeline in `1-network-centrality/processed-data`: `{output_prefix}-CentralityCosDist-input.csv` and `{output_prefix}-{organism_label}-nodes-centrality.csv`. 

The values of {output_prefix} and {organism_label} are taken from the configuration file at runtime. 

`{output_prefix}-CentralityCosDist-input.csv` is an intermediate file written in the specific format required for Centrality Cosine Distance calculations. 

`{output_prefix}-{organism_label}-step1.csv` is the output file from this step containing the annotated node network. Each row corresponds to a node, with columns for centrality metrics added for each node. 

Centrality Cosine Distance calculations are carried out using code from the [Mukhtar Lab](https://github.com/nilesh-iiita/CentralityCosDist). This repository is distributed in `1-network-centrality/python-scripts/CentralityCosDist`. 

The directory `1-network-centrality/test` contains input data used to perform a test at runtime that the CentralityCosDist code is functioning as expected. 

Note: In the current version of the repository, CentralityCosDist calculations are run but not included in the data product. They can be added back in by commenting out the line

```python
nodes_df = nodes_df.drop(columns=["CentralityCosDist_rank", "CentralityCosDist_similarity_score"])
```
 in `python-scripts/network-centrality.py`
