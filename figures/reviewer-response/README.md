The program `define-missing-nodes.ipynb` is used to define the set of nodes that do not occur in either the AE-MS or Y2H networks, `processed-data/missing-nodes.csv` which is used as an input by the configuration file `../../config-files/missing.config` to annotate all nodes in neither the AE-MS or Y2H networks. 

The program `make-missing-data-files.ipynb` uses the annotated AE-MS, Y2H, and "missing" pkl files to generate to pd.DataFrame objects and save them as pickle and csv:

`processed-data/missing-from-y2h-s288c.*`
`processed-data/missing-from-ae-ms-s288c.*`

These files can then be used straightforwardly to compare in-network and missing values. 
