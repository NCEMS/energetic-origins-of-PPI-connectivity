### Annotate network with centrality metrics

Computes the following centrality metrics:

* Degree centrality
* Betweenness centrality
* Eigenvector centrality
* Closeness centrality
* Load centrality
* Page rank
* Information centrality (run on the largest connected subgraph)
* CentralityCosDist (using each of the above as part of the vector; default is to use all nodes as seed nodes)

Weighted k-shell decomposition is pre-computed in Cytoscape using The Yeast Interactome as the input; this information is stored in ../0-download-inputs/data-files/The_Yeast_Interactome_nodes.csv
