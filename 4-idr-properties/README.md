### Predict IDRs and annotate them with sequence- and ensemble-based properties

* Run `metapredict` on the `trimmed_sequence` of each protein and extract IDR sequences
* For each IDR sequence, run Albatross to extract properties of the conformation ensemble
* For each IDR sequence, run CIDER to extract sequence-based properties

To run this pipeline step in isolation, run the command `snakemake --use-conda` from the `4-idr-properties` directory.
