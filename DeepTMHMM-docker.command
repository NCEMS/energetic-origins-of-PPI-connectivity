mkdir deeptmhmm-results

FASTA=$(realpath data-files/partial-orf_trans.fasta)
TMPDIR=$(realpath tmp)
OUTDIR=$(realpath deeptmhmm-results)

touch $OUTDIR/predicted_topologies.3line
touch $OUTDIR/TMRs.gff3
touch $OUTDIR/deeptmhmm_results.md
touch $OUTDIR/plot.png

docker run  \
    -v $FASTA:/openprotein/prot_seqs.fasta \
    -v $TMPDIR/embeddings:/openprotein/embeddings \
    -v $TMPDIR/probabilities:/openprotein/probabilities \
    -v $OUTDIR/predicted_topologies.3line:/openprotein/predicted_topologies.3line \
    -v $OUTDIR/TMRs.gff3:/openprotein/TMRs.gff3 \
    -v $OUTDIR/deeptmhmm_results.md:/deeptmhmm_results.md \
    -v $OUTDIR/plot.png:/openprotein/plot.png \
    biswasaneel/deeptmhmm:latest \
    python3 predict.py --fasta /openprotein/prot_seqs.fasta
