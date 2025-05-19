from pathlib import Path

LOCAL_CONFIG   = config["translation_speed"]

DATA_DIR       = Path(LOCAL_CONFIG["input_dir"])
WORK_DIR       = Path(workflow.basedir)
OUTPUT_PREFIX  = config["output_prefix"]
OUTPUT_DIR     = WORK_DIR / LOCAL_CONFIG["output_dir"]
ORGANISM_TAG   = config["organism_label"]
RAW_DATA_DIR   = Path(LOCAL_CONFIG["raw_data_dir"])
RIBO_FILES     = [f"{RAW_DATA_DIR}/{fname}" for fname in LOCAL_CONFIG["ribo_seq_data"]]

rule all:
    input:
        f"{OUTPUT_DIR}/{OUTPUT_PREFIX}-{ORGANISM_TAG}-nodes-centrality-seqs-DeepTMHMM-SignalP-UniProt-IDRs-albatross-cider-GhoshDill-Cagiada-Rosetta-FoldX-halflife-expr-speed.pkl"

rule prepare_reads:
    input:
        data = RIBO_FILES
    conda:
        "env/add-translation-speed.yml"
    params:
        output_prefix = OUTPUT_PREFIX,
        output_dir    = str(OUTPUT_DIR),
        organism_tag  = ORGANISM_TAG,
        ribo_files    = " ".join(RIBO_FILES)
    output:
        f"{OUTPUT_DIR}/{OUTPUT_PREFIX}-{ORGANISM_TAG}-pooled-ribo-seq-data.pkl"
    shell:
        f"""
        python {WORK_DIR}/python-scripts/prepare-reads.py \
        --output_dir {{params.output_dir}} \
        --output_prefix {{params.output_prefix}} \
        --organism_tag {{params.organism_tag}} \
        --ribo_seq_data_files {{params.ribo_files}}
        """

rule add_translation_speed:
    input:
        nodes = f"{DATA_DIR}/{OUTPUT_PREFIX}-{ORGANISM_TAG}-nodes-centrality-seqs-DeepTMHMM-SignalP-UniProt-IDRs-albatross-cider-GhoshDill-Cagiada-Rosetta-FoldX-halflife-expr.pkl",
        ta_data = f"{OUTPUT_DIR}/{OUTPUT_PREFIX}-{ORGANISM_TAG}-pooled-ribo-seq-data.pkl"
    conda:
        "env/add-translation-speed.yml"
    params:
        output_prefix  = OUTPUT_PREFIX,
        output_dir     = str(OUTPUT_DIR),
        organism_tag   = ORGANISM_TAG
    output:
        f"{OUTPUT_DIR}/{OUTPUT_PREFIX}-{ORGANISM_TAG}-nodes-centrality-seqs-DeepTMHMM-SignalP-UniProt-IDRs-albatross-cider-GhoshDill-Cagiada-Rosetta-FoldX-halflife-expr-speed.pkl"
    shell:
        f"""
        python {WORK_DIR}/python-scripts/add-translation-speed.py \
        --nodes            {{input.nodes}} \
        --output_dir       {{params.output_dir}} \
        --output_prefix    {{params.output_prefix}} \
        --organism_tag     {{params.organism_tag}} \
        --trans_speed_file {{input.ta_data}}
        """
