from pathlib import Path
import os
from glob import glob

# === CONFIG ===
LOCAL_CONFIG      = config["dG_predictions"]
WORK_DIR          = Path(workflow.basedir)
DATA_DIR1         = Path(LOCAL_CONFIG["input_dir"])
OUTPUT_PREFIX     = config["output_prefix"]
OUTPUT_DIR        = WORK_DIR / LOCAL_CONFIG["output_dir"]
ORGANISM_TAG      = config["organism_label"]
TEMPERATURE       = LOCAL_CONFIG["temperature"]
RUN_CAGIADA       = LOCAL_CONFIG["run_cagiada"]
STRUCTURE_DIR     = Path(LOCAL_CONFIG["structure_dir"])
ESM_MODEL         = Path(LOCAL_CONFIG["ESM_model"])
SEQ_COLUMN        = LOCAL_CONFIG["seq_column"]
TEST_DIR          = WORK_DIR / LOCAL_CONFIG["test_dir"]
AF2_DATABASES     = LOCAL_CONFIG["AF2_databases"]
MAX_TEMPLATE_DATE = LOCAL_CONFIG["max_template_date"]

# === GPU assignment ===
def get_gpu_for_node(node):
    return int(hash(node)) % 2

# === Checkpoint-based target discovery ===
def discover_alphafold_targets(checkpoint_output_dir):
    base_dir = os.path.join(checkpoint_output_dir, "alphafold")
    targets = []
    for path in sorted(glob(f"{base_dir}/*/*.fasta")):
        node = os.path.basename(path).replace(".fasta", "")
        targets.append(node)
    return targets

def get_af2_outputs(wildcards):
    checkpoint_output = checkpoints.prepare_structures_checkpoint.get()
    fasta_dir = checkpoint_output.output.fasta_dir
    targets = discover_alphafold_targets(fasta_dir)
    return expand(
        f"{config['dG_predictions']['output_dir']}/alphafold/{{node}}/result_model_1_pred_0.pdb",
        node=targets
    )

# === RULES ===

rule all:
    input:
        get_af2_outputs,
        str(OUTPUT_DIR / f"{OUTPUT_PREFIX}-{ORGANISM_TAG}-nodes-centrality-seqs-DeepTMHMM-SignalP-UniProt-IDRs-albatross-cider-GhoshDill-Cagiada.pkl")


checkpoint prepare_structures_checkpoint:
    input:
        nodes = str(DATA_DIR1 / f"{OUTPUT_PREFIX}-{ORGANISM_TAG}-nodes-centrality-seqs-DeepTMHMM-SignalP-UniProt-IDRs-albatross-cider.pkl")
    output:
        struct_pkl = str(OUTPUT_DIR / f"{OUTPUT_PREFIX}-{ORGANISM_TAG}-temp.pkl"),
        fasta_dir  = directory(f"{config['dG_predictions']['output_dir']}/alphafold")
    conda:
        "env/prepare-structures.yml"
    params:
        input_dir     = str(STRUCTURE_DIR),
        output_dir    = str(OUTPUT_DIR),
        organism_tag  = ORGANISM_TAG,
        output_prefix = OUTPUT_PREFIX
    shell:
        f"""
        python {WORK_DIR}/python-scripts/prepare-structures.py \
        --nodes         {{input.nodes}} \
        --input_dir     {{params.input_dir}} \
        --output_dir    {{params.output_dir}} \
        --output_prefix {{params.output_prefix}} \
        --organism_tag  {{params.organism_tag}}
        """


rule run_alphafold2:
    input:
        fasta = lambda wildcards: f"{config['dG_predictions']['output_dir']}/alphafold/{wildcards.node}/{wildcards.node}.fasta",
        struct_pkl = str(OUTPUT_DIR / f"{OUTPUT_PREFIX}-{ORGANISM_TAG}-temp.pkl")
    output:
        pdb = f"{config['dG_predictions']['output_dir']}/alphafold/{{node}}/result_model_1_pred_0.pdb"
    params:
        data_dir = AF2_DATABASES,
        max_template_date = MAX_TEMPLATE_DATE,
        gpu_id = lambda wildcards: get_gpu_for_node(wildcards.node)
    shell:
        """
        docker run --rm --gpus all \
          -v {os.path.dirname(input.fasta)}:/mnt/fasta_path_0:ro \
          -v {params.data_dir}:/mnt/data_dir:ro \
          -v {os.path.dirname(output.pdb)}:/mnt/output \
          -e NVIDIA_VISIBLE_DEVICES={params.gpu_id} \
          -e TF_FORCE_UNIFIED_MEMORY=1 \
          -e XLA_PYTHON_CLIENT_MEM_FRACTION=0.8 \
          --user $(id -u):$(id -g) \
          alphafold \
          --fasta_paths=/mnt/fasta_path_0/{wildcards.node}.fasta \
          --output_dir=/mnt/output \
          --data_dir=/mnt/data_dir \
          --max_template_date={params.max_template_date} \
          --model_preset=monomer \
          --db_preset=full_dbs \
          --uniref90_database_path=/mnt/data_dir/uniref90/uniref90.fasta \
          --mgnify_database_path=/mnt/data_dir/mgnify/mgy_clusters_2022_05.fa \
          --template_mmcif_dir=/mnt/data_dir/pdb_mmcif/mmcif_files \
          --obsolete_pdbs_path=/mnt/data_dir/pdb_mmcif/obsolete.dat \
          --pdb70_database_path=/mnt/data_dir/pdb70/pdb70 \
          --uniref30_database_path=/mnt/data_dir/uniref30/UniRef30_2021_03 \
          --bfd_database_path=/mnt/data_dir/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt \
          --use_gpu_relax=False
        """

rule ghosh_dill_stability:
    input:
        nodes = str(OUTPUT_DIR / f"{OUTPUT_PREFIX}-{ORGANISM_TAG}-temp.pkl")
    conda:
        "env/ghosh-dill-stability.yml"
    params:
        output_prefix  = OUTPUT_PREFIX,
        output_dir     = str(OUTPUT_DIR),
        organism_tag   = ORGANISM_TAG,
        temperature    = TEMPERATURE,
        seq_column     = SEQ_COLUMN
    output:
        nodes          = str(OUTPUT_DIR / f"{OUTPUT_PREFIX}-{ORGANISM_TAG}-nodes-centrality-seqs-DeepTMHMM-SignalP-UniProt-IDRs-albatross-cider-GhoshDill.pkl")
    shell:
        f"""
        python {WORK_DIR}/python-scripts/ghosh-dill-stability.py \
        --nodes         {{input.nodes}} \
        --output_dir    {{params.output_dir}} \
        --output_prefix {{params.output_prefix}} \
        --organism_tag  {{params.organism_tag}} \
        --temperature   {{params.temperature}} \
        --seq_column_to_use {{params.seq_column}}
        """

rule cagiada_stability:
    input:
        nodes = str(OUTPUT_DIR / f"{OUTPUT_PREFIX}-{ORGANISM_TAG}-nodes-centrality-seqs-DeepTMHMM-SignalP-UniProt-IDRs-albatross-cider-GhoshDill.pkl")
    conda:
        "env/cagiada-stability.yml"
    params:
        output_prefix  = OUTPUT_PREFIX,
        output_dir     = str(OUTPUT_DIR),
        organism_tag   = ORGANISM_TAG,
        structure_dir  = STRUCTURE_DIR,
        ESM_model      = ESM_MODEL,
        test_dir       = TEST_DIR,
        run_cagiada    = RUN_CAGIADA
    output:
        str(OUTPUT_DIR / f"{OUTPUT_PREFIX}-{ORGANISM_TAG}-nodes-centrality-seqs-DeepTMHMM-SignalP-UniProt-IDRs-albatross-cider-GhoshDill-Cagiada.pkl")
    shell:
        f"""
        python {WORK_DIR}/python-scripts/cagiada-stability.py \
        --nodes          {{input.nodes}} \
        --output_dir     {{params.output_dir}} \
        --output_prefix  {{params.output_prefix}} \
        --organism_tag   {{params.organism_tag}} \
        --structure_dir  {{params.structure_dir}} \
        --ESM_model      {{params.ESM_model}} \
        --test_dir       {{params.test_dir}} \
        --run_cagiada    {{params.run_cagiada}} \
        > {WORK_DIR}/cagiada.log
        """
