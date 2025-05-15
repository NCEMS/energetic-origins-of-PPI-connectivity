from pathlib import Path

LOCAL_CONFIG  = config["uniprot_annotation"]

WORK_DIR      = Path(workflow.basedir)
DATA_DIR1     = Path(LOCAL_CONFIG["input_dir"])
DATA_DIR2     = Path(config["data_dir"])
OUTPUT_PREFIX = config["output_prefix"]
OUTPUT_DIR    = WORK_DIR / LOCAL_CONFIG["output_dir"]
ORGANISM      = LOCAL_CONFIG["organism"]
ORGANISM_TAG  = config["organism_label"]

rule all:
    input:
        str(OUTPUT_DIR / f"{ORGANISM_TAG}-uniprot_sprot.csv"),
        str(OUTPUT_DIR / f"{OUTPUT_PREFIX}-{ORGANISM_TAG}-nodes-centrality-seqs-DeepTMHMM-SignalP-UniProt.csv")

rule parse_uniprot:
    input:
        uniprot = str(DATA_DIR2 / "uniprot_sprot.dat")
    conda:
        "env/parse-uniprot.yml"
    params:
        output_dir   = str(OUTPUT_DIR),
        organism     = ORGANISM,
        output_file  = str(OUTPUT_DIR / f"{ORGANISM_TAG}-uniprot_sprot.csv")
    output:
        str(OUTPUT_DIR / f"{ORGANISM_TAG}-uniprot_sprot.csv")
    shell:
        f"""
        python {WORK_DIR}/python-scripts/parse-uniprot.py \
        --input_file    {{input.uniprot}} \
        --output_file   {{params.output_file}} \
        --organism      "{{params.organism}}"
        """

rule add_uniprot_info:
    input:
        nodes = str(DATA_DIR1 / f"{OUTPUT_PREFIX}-{ORGANISM_TAG}-nodes-centrality-seqs-DeepTMHMM-SignalP.csv"),
        unipr = str(OUTPUT_DIR / f"{ORGANISM_TAG}-uniprot_sprot.csv")
    conda:
        "env/add-uniprot-info.yml"
    params:
        output_prefix = OUTPUT_PREFIX,
        output_dir    = str(OUTPUT_DIR),
        organism_tag  = ORGANISM_TAG
    output:
        str(OUTPUT_DIR / f"{OUTPUT_PREFIX}-{ORGANISM_TAG}-nodes-centrality-seqs-DeepTMHMM-SignalP-UniProt.csv")
    shell:
        f"""
        python {WORK_DIR}/python-scripts/add-uniprot-info.py \
        --nodes         {{input.nodes}} \
        --uniprot       {{input.unipr}} \
        --output_prefix {{params.output_prefix}} \
        --output_dir    {{params.output_dir}} \
        --organism_tag  {{params.organism_tag}}
        """
