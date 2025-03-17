rule download_inputs:
	output:
		"data-files/.download_complete"
	shell:
		"""
		bash bash_scripts/download_inputs.sh
		touch {output}
		"""
