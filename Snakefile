# Define global variables for generalization
OUTPUT_PREFIX = "0_"
DATA_DIR = "data-files"
PROCESSED_DIR = "processed-data"

rule all:
    input:
        f"{DATA_DIR}/.all_fasta_created",
        f"{PROCESSED_DIR}/{OUTPUT_PREFIX}network_nodes_with_annotation.csv",
        f"{PROCESSED_DIR}/{OUTPUT_PREFIX}network_nodes_with_stability.csv"

rule download_inputs:
    params:
        fDir = DATA_DIR
    output:
        f"{DATA_DIR}/orf_trans.fasta",
        f"{DATA_DIR}/YEAST_559292_idmapping.dat"
    shell:
        """
        bash bash-scripts/download-inputs.sh {params.fDir}
        """

rule create_fasta:
    params:
        fDir = DATA_DIR
    output:
        f"{DATA_DIR}/.all_fasta_created"
    shell:
        """
        bash bash-scripts/extract-seqres.sh {params.fDir}
        """

rule network_analysis:
    input:
        fasta = f"{DATA_DIR}/orf_trans.fasta",
        edges = f"{DATA_DIR}/The_Yeast_Interactome_edges.csv",
        preds = "DeepTMHMM-runs/s288c-results/all-predictions-s288c.3line"
    params:
        output_prefix = OUTPUT_PREFIX,
        output_dir = PROCESSED_DIR
    output:
        f"{PROCESSED_DIR}/{OUTPUT_PREFIX}network_nodes_with_annotation.csv"
    shell:
        """
        conda run -n network-analysis python python-scripts/network-analysis.py \
        --edges {input.edges} --fasta {input.fasta} --output_prefix {params.output_prefix} \
        --output_dir {params.output_dir} --seq_preds {input.preds}
        """

rule cagiada_stability:
    input:
        input_node_file = f"{PROCESSED_DIR}/{OUTPUT_PREFIX}network_nodes_with_annotation.csv"
    params:
        output_dir = PROCESSED_DIR
    output:
        f"{PROCESSED_DIR}/{OUTPUT_PREFIX}network_nodes_with_stability.csv"
    shell:
        """
        conda run -n cagiada-stability python python-scripts/cagiada-stability.py \
        --input_node_file {input.input_node_file} --output_dir {params.output_dir}
        """
