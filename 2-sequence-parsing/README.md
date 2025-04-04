### Add sequence information and DeepTMHMM annotations to nodes

Add sequence information from `orf_trans.fasta` and DeepTMHMM annotations for each sequence. 

Execute this portion of the pipeline by running the command

`snakemake --use-conda`

Note: signal sequences are currently cleaved off of the sequence based on DeepTMHMM predictions. This part of the code will be updated to include cleavage sites based on SignalP results ASAP. 
