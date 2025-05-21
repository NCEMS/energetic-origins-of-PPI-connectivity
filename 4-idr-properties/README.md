### Predict IDRs and annotate them with sequence- and ensemble-based properties

* Run metapredict v3.0 on the `trimmed_sequence` of each protein and extract IDR sequences
* For each IDR sequence, run ALBATROSS to extract properties of the conformation ensemble
* For each IDR sequence, run CIDER to extract sequence-based properties (both ALBATROSS and CIDER can be run from SPARROW)

*Note*: The IDR section of the configuration files contains the parameters `min_idr_length`, which determines how large a contiguous set of residues needs to be before it "counts" as an IDR, and `disorder_resi_cut`, which determines the metapredict disorder threshold above which a residue is considered to be disordered. Do not change `disorder_resi_cut` without careful consideration, as this is the default within metapredict v3.0
