### `flatten`: Finally, produce the analysis-ready data product

Corresponding configuration file section: `Flatten nodes`

* This pipeline step performs some simple cleanup of your final data product, like removing duplicate column names
* As its name implies, this pipeline step also produces a version of the database which is "flattened" over IDRs, produced a version of the annotated network with one row per IDR rather than one row per protein
* If you are running with data other than The Yeast Interactome and associated annotations, we recommend writing a new `python-scripts/flatten.py` program to manage final cleanup of your data product
