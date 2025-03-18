rule download_inputs:

	output:
		"data-files/.download_complete"

	shell:
		"""
		bash bash_scripts/download_inputs.sh
		touch {output}
		"""

rule network_analysis

	input:
		fasta="data-files/orf_trans.fasta"
		edges="data-files/The_Yeast_Interactome_edges.csv"
		output_prefix="0_"
		output_dir="processed-data"

	output:
		"processed-data/0_network_nodes_with_annotation.csv

	shell:
		"""
		conda run -n network-analysis python python-scripts/network-analysis.py --edges {edges} --fasta {fasta} --output_prefix {output_prefix} --output_dir {output_dir}
		"""

rule cagiada_stability

	input:
		

	output:

	shell:

		"""
		conda run -n cagiada-stability python python-scripts/cagiada-stability.py 
		"""
