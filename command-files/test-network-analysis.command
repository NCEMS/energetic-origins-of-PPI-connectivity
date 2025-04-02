python python-scripts/network-analysis.py \
--edges test/network/data-files/edges.csv \
--fasta data-files/orf_trans.fasta \
--structure_dir data-files \
--output_prefix test_ \
--output_dir test/network/processed-data \
--seq_preds DeepTMHMM-runs/s288c-results/all-predictions-s288c.3line \
--ID_mappings data-files/YEAST_559292_idmapping.dat \
--uniprot_data processed-data/uniprot_sprot-s288c.csv
