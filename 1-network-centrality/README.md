### Annotate network with centrality metrics

Computes the following centrality metrics:

* Degree centrality
* Betweenness centrality
* Eigenvector centrality
* Closeness centrality
* Load centrality
* Page rank
* k-shell (core number) - *To be changed to weighted k-shell*
* CentralityCosDist (using each of the above as part of the vector; default is to use all nodes as seed nodes)

To run this pipeline step in isolation, run the command `snakemake --use-conda` from the `1-network-centrality` directory.

`snakemake --use-conda`
