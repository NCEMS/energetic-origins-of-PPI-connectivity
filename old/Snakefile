# global variables
OUTPUT_PREFIX = "0_"
DATA_DIR      = "/home/jovyan/data-store/home/shared/NCEMS/working-groups/energetic-origins/data-files"
PROCESSED_DIR = "/home/jovyan/data-store/home/shared/NCEMS/working-groups/energetic-origins/processed-data"

# rule that sets overall outputs required by this pipeline
rule all:
    input:
        f"{DATA_DIR}/uniprot_sprot.dat",
        f"{DATA_DIR}/.all_fasta_created",
        f"{PROCESSED_DIR}/{OUTPUT_PREFIX}network_nodes_with_annotation.csv",
        f"{PROCESSED_DIR}/{OUTPUT_PREFIX}network_nodes_with_annotation_and_stability.csv",
        f"{PROCESSED_DIR}/{OUTPUT_PREFIX}annotated_network_summary.csv"

# download all inputs required (other than The Yeast Interactome files)
rule download_inputs:
    params:
        fDir = DATA_DIR
    output:
        f"{DATA_DIR}/orf_trans.fasta",
        f"{DATA_DIR}/YEAST_559292_idmapping.dat",
        f"{DATA_DIR}/SGD_features.tab",
        f"{DATA_DIR}/uniprot_sprot.dat"
    shell:
        """
        bash bash-scripts/download-inputs.sh {params.fDir}
        """

# extract SEQRES records from AF2 PDB files
rule create_fasta:
    params:
        fDir = DATA_DIR
    output:
        f"{DATA_DIR}/.all_fasta_created"
    shell:
        """
        bash bash-scripts/extract-seqres.sh {params.fDir}
        """

# parse information from UniProt
rule parse_uniprot:
    params:
        uniprot_db = f"{DATA_DIR}/uniprot_sprot.dat"
    output:
        f"{PROCESSED_DIR}/uniprot_sprot-s288c.csv"
    shell:
        """
        conda run -n network-analysis python python-scripts/parse-uniprot.py \
        --input_file params.uniprot_db
        --output_file f"{PROCESSED_DIR}/uniprot_sprot-s288c.csv"
        --organism "Saccharomyces cerevisiae (strain ATCC 204508 / S288c) (Baker's yeast)."
        """

# annotate network with centrality metrics, IDR information (metapredict v3)
# and use DeepTMHMM to annotate membrane proteins, signal peptides, etc. 
rule network_analysis:
    input:
        fasta = f"{DATA_DIR}/orf_trans.fasta",
        edges = f"{DATA_DIR}/The_Yeast_Interactome_edges.csv",
        preds = "DeepTMHMM-runs/s288c-results/all-predictions-s288c.3line",
        mappi = f"{DATA_DIR}/YEAST_559292_idmapping.dat",
        unipr = f"{PROCESSED_DIR}/uniprot_sprot-s288c.csv"
    params:
        output_prefix = OUTPUT_PREFIX,
        output_dir    = PROCESSED_DIR,
        structure_dir = DATA_DIR
    output:
        f"{PROCESSED_DIR}/{OUTPUT_PREFIX}network_nodes_with_annotation.csv"
    shell:
        """
        conda run -n network-analysis python python-scripts/network-analysis.py \
        --edges         {input.edges} \
        --fasta         {input.fasta} \
        --output_prefix {params.output_prefix} \
        --output_dir    {params.output_dir} \
        --ID_mappings   {input.mappi} \
        --structure_dir {params.structure_dir} \
        --seq_preds     {input.preds} \
        --uniprot_data  {input.unipr}
        """

# add stability predictions to nodes
rule compute_dG:
    input:
        input_node_file = f"{PROCESSED_DIR}/{OUTPUT_PREFIX}network_nodes_with_annotation.csv"
    params:
        output_prefix = OUTPUT_PREFIX,
        output_dir = PROCESSED_DIR
    output:
        f"{PROCESSED_DIR}/{OUTPUT_PREFIX}network_nodes_with_annotation_and_stability.csv"
    shell:
        """
        conda run -n cagiada-stability python python-scripts/compute-dG.py \
        --input_node_file {input.input_node_file} \
        --output_dir      {params.output_dir} \
        --output_prefix   {params.output_prefix} \
        --do_cagiada      False \
        --temperature     303.15
        """

# profile the output DataFrame
rule profile_output:
    input:
        input_annotated_network = f"{PROCESSED_DIR}/{OUTPUT_PREFIX}network_nodes_with_annotation_and_stability.csv"
    params:
        output_prefix = OUTPUT_PREFIX,
        output_dir    = PROCESSED_DIR
    output:
        f"{PROCESSED_DIR}/{OUTPUT_PREFIX}annotated_network_summary.csv"
    shell:
        """
        conda run -n network-analysis python python-scripts/profile-output.py \
        --input_annotated_network {input.input_annotated_network} \
        --output_prefix           {params.output_prefix} \
        --output_dir              {params.output_dir}
        """
