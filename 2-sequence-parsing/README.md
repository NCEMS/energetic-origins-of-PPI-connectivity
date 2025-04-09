### Add sequence information and DeepTMHMM annotations to nodes

Add sequence information from `orf_trans.fasta` and DeepTMHMM annotations for each sequence. 

To run this pipeline step in isolation, run the command `snakemake --use-conda` from the `2-sequence-parsing` directory.

*N.B.*: signal sequences are currently cleaved based on DeepTMHMM predictions. This part of the code will be updated to include cleavage sites based on SignalP results ASAP. 
