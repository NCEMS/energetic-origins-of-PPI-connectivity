### `4-idr-properties`: Predict IDRs and annotate them with sequence- and ensemble-based properties

Corresponding configuration file section: `IDR IDENTIFICATION & PROPERTY PREDICTION`

Anticipated execution time: ~10 min

To run this pipeline step in isolation, run the command:

```bash
snakemake -c all --use-conda --conda-frontend conda --snakefile 4-idr-properties/Snakefile --configfile config-files/s288c.config
```

* Run metapredict v3.0 on the `signalP_trimmed_sequence` of each protein to extract IDR sequences
* Classify proteins as disordered or not (binary) using DisProt-based threshold(s)
* For each IDR sequence, run ALBATROSS to extract properties of the conformation ensemble
* For each IDR sequence, run CIDER to extract sequence-based properties (both ALBATROSS and CIDER can be run from SPARROW)

*Note*: The IDR section of the configuration files contains the parameter `min_idr_length`, which determines how large a contiguous set of residues needs to be before it "counts" as a discrete IDR. The number of IDRs within a protein will change depending on this threshold. 
