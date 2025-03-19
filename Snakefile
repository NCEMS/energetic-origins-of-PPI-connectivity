rule all:
	input:
		"processed-data/0_network_nodes_with_annotation.csv",
		"processed-data/.dG_done"

rule download_inputs:
	output:
		"data-files/orf_trans.fasta",
		"data-files/YEAST_559292_idmapping.dat"
	shell:
		"""
		bash bash-scripts/download-inputs.sh
		"""

rule network_analysis:
	input:
		fasta           = "data-files/orf_trans.fasta",
		edges           = "data-files/The_Yeast_Interactome_edges.csv"
	params:
		output_prefix   = "0_",
		output_dir      = "processed-data"
	output:
		"processed-data/0_network_nodes_with_annotation.csv"
	shell:
		"""
		conda run -n network-analysis python python-scripts/network-analysis.py \
		--edges {input.edges} --fasta {input.fasta} --output_prefix {params.output_prefix} \
		--output_dir {params.output_dir}
		"""

rule cagiada_stability:
	input:
		input_node_file = "processed-data/0_network_nodes_with_annotation.csv"
	params:
		output_dir      = "processed-data"
	output:
		"processed-data/.dG_done"
	shell:
		"""
		conda run -n cagiada-stability python python-scripts/cagiada-stability.py \
		--input_node_file {input.input_node_file} --output_dir {params.output_dir}
		touch {output}
		"""
