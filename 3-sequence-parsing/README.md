### `2-sequence-parsing`: Add sequence information and DeepTMHMM/SignalP annotations to nodes

Corresponding configuration file section: `SEQUENCE PARSING`

Critical notes for understanding the A. thaliana processing pipeline:

Node/edge data do not include isoform IDs, but the TAIR12 proteome release used to source sequence information does contain isoform IDs. We always choose isoform 1, suffix ".1", from the TAIR12 proteome sequences. 

Below is the current stdout from this pipeline step:

```text
Warning: some nodes had TAIR12 sequence records but no matching .1 isoform.
Number of affected nodes: 1
Examples:
AT1G06515: expected AT1G06515.1; found AT1G06515.2
Input nodes: 11,695
Nodes with TAIR12 .1 sequence: 11,550
Nodes without TAIR12 .1 sequence: 145
AGIs with at least one UniProtKB-AC candidate: 28,613
AGIs with multiple UniProtKB-AC candidates: 10,362
UniProt FASTA accessions parsed: 27,496

UniProt sequence-mapping status counts:
uniprot_match_status
exact_sequence_match                      9058
best_nonexact_sequence_match              1570
no_candidate_sequence_in_uniprot_fasta     906
no_tair12_sequence                         119
no_uniprot_candidates                       42

UniProt exact sequence-match counts:
uniprot_exact_sequence_match
True     9058
False    1570
<NA>     1067

Wrote sequence-annotated nodes to: /home/dan182/energetic-origins/arabidopsis-branch/energetic-origins-of-PPI-connectivity/3-sequence-parsing/processed-data/union-athaliana-seqs-step3.csv
```

Anticipated execution time: ~5 min

To run this pipeline step in isolation, run the command:

```bash
snakemake -c all --use-conda --conda-frontend conda --snakefile 2-sequence-parsing/Snakefile --configfile config-files/union-athaliana.config
```

### LEGACY S288C README.md


```bash
snakemake -c all --use-conda --conda-frontend conda --snakefile 2-sequence-parsing/Snakefile --configfile config-files/s288c.config
```

This pipeline will:

(1) Add sequence information for genes as possible (sourced from SGD; see `../0-download-inputs/README.md` for details)
(2) Add DeepTMHMM predictions for each sequence (DeepTMHMM currently must be run outside of the pipeline and the path to the results included in the .config)
(3) Predict signal peptide cleavage sites with SignalP
(4) Add updated sequence information in "*_trimmed_sequence" column(s) of the annotated node network

* SignalP, when able to access a GPU, will require ~15 min. All other rules in this step require <2 min. 

* To run SignalP after downloading and unpacking the directory as described in `../README.md`, follow these steps:

* You *must* copy the model weights from the model folder into the signalP folder:

`cp 2-sequence-parsing/python-scripts/signalp6_fast/signalp-6-package/models/distilled_model_signalp6.pt 2-sequence-parsing/python-scripts/signalp6_fast/signalp-6-package/signalp/model_weights/`

* You *must* edit the environment file `env/signalp_fast_local.yml` to insert the absolute path to the directory `signalp-6-package`. In the version distributed with this repository, the path is `/home/dan182/energetic-origins/energetic-origins-of-PPI-connectivity/2-sequence-parsing/python-scripts/signalp6_fast/signalp-6-package`

DeepTMHMM notes:

* The directory `DeepTMHMM-runs` contains the pre-computed results for the S288C yeast open reading frames.
* The file `DeepTMHMM-runs/command` provides the command used to generate the results in `DeepTMHMM-runs/s288c-results/all-predictions-s288c.3line`
* To rerun this section, install DockerHub and the image `dtu/deeptmhmm:1.0.24`. Start the container and then run the command in `DeepTMHMM-runs/command`
* Processing all sequences will require 24-36 hours without a GPU. 
