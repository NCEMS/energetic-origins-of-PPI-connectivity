### `0-download-inputs`: Download input files used by downstream pipeline steps

This section of the pipeline handles the download of required files to `data-files/`:

The first rule of the Snakemake pipeline will download the files.

After downloading the files, the second rule will extract SEQRES records from each AF2 PDB file and convert them to fasta format

This step will require 60-90 min depending on connection speeds and the write speed of your file system. 

*N.B.*: Once all files are downloaded and unpacked, `0-download-inputs/data-files` will contain ~112 GB of data. 

#### Summary of data files

Some files are distributed with the GitHub repository while others are downloaded when running the Snakemake pipeline. 

After running this step, you should have the following file counts:

| File extension | Number of files |
|:--------------:|----------------:|
| pdb            | 6,168           |
| fasta          | 6,040           |
| csv            |    20           |
| xlsx           |     9           |
| txt            |     5           |
| tsv            |     4           |
| dat            |     3           |
| obo            |     1           |
| pt             |     1           |
| xml            |     1           |
| cys            |     1           |
| TOTAL FILES    |12,253           |


The following files are provided with the GitHub repository:

|Used in step| Filename | Description | Source |
|:----------:|:---------|:------------|:-------|

|  1 | The_Yeast_Interactome.cys | The Yeast Interactome Cytoscape session file | Downloaded from https://www.yeast-interactome.org/ | 
|  1 | The_Yeast_Interactome_nodes.csv | List of nodes (proteins) | Extracted from .cys file with Cytoscape |
|  1 | The_Yeast_Interactome_edges.csv | List of edges (interactions) between nodes | Extracted from .cys file with Cytoscape |
|  4 | DisProt-release_2024_12-with_ambiguous_evidences.tsv | DisProt database version 2024_12 | Downloaded from https://disprot.org/download |
|  4 | DisProt-release_2024_12-with_ambiguous_evidences-cleaned.tsv | DisProt database version 2024_12 | Cleaned version of DisProt with whitespace cleaned for reading |
|  8 | 1-s2.0-S2211124714009346-mmc2.xlsx | Protein half-lives in S. cerevisiae | 10.1016/j.celrep.2014.10.065, Table S1 |
|  8 | 1-s2.0-S2211124714009346-mmc2.csv | Protein half-lives in S. cerevisiae | 10.1016/j.celrep.2014.10.065, Table S1 (reformatted for analysis) |
|  9 | 1-s2.0-S240547121730546X-mmc5.xlsx | Protein expression in S. cerevisiae | 10.1016/j.cels.2017.12.004, Table S4 |
|  9 | 1-s2.0-S240547121730546X-mmc5.csv | Protein expression in S. cerevisiae | 10.1016/j.cels.2017.12.004, Table S4 (reformatted for analysis)|
| 10 | skr_weinberg_genesTE.csv | Translation efficiency information from scikit-ribo for S. cerevisiae | Downloaded from https://github.com/schatzlab/scikit-ribo_manuscript/blob/master/Data/skr_weinberg_genesTE.csv |
| 12 | 20210721_YeastRefold_ProteinSummary_1min_meta.xlsx | LiP-MS data; 1 min after refolding initiated, sample 1 | Fried Lab, unpublished |
| 12 | 20210721_YeastRefold_ProteinSummary_1min_meta.csv | LiP-MS data; 1 min after refolding initiated, sample 1 (reformatted from matched .xlsx) | Fried Lab, unpublished |
| 12 | 20220505_YeastRefolding_ProteinSummary_1min_meta.xlsx | LiP-MS data; 1 min after refolding initiated, sample 2 | Fried Lab, unpublished |
| 12 | 20220505_YeastRefolding_ProteinSummary_1min_meta.csv | LiP-MS data; 1 min after refolding initiated, sample 2 (reformatted from matched .xlsx) | Fried Lab, unpublished |
| 12 | 20220505_YeastRefolding_ProteinSummary_1min_duplicate_meta.xlsx | LiP-MS data; 1 min after refolding initiated, sample 3 | Fried Lab, unpublished |
| 12 | 20220505_YeastRefolding_ProteinSummary_1min_duplicate_meta.csv | LiP-MS data; 1 min after refolding initiated, sample 3 (reformatted from matched .xlsx) | Fried Lab, unpublished |
| 12 | 20210721_YeastRefold_ProteinSummary_5min_meta.xlsx | LiP-MS data; 5 min after refolding initiated, sample 1 | Fried Lab, unpublished |
| 12 | 20210721_YeastRefold_ProteinSummary_5min_meta.csv | LiP-MS data; 5 min after refolding initiated, sample 1 (reformatted from matched .xlsx) | Fried Lab, unpublished |
| 12 | 20220505_YeastRefolding_ProteinSummary_5min_meta.xlsx | LiP-MS data; 5 min after refolding initiated, sample 2 | Fried Lab, unpublished |
| 12 | 20220505_YeastRefolding_ProteinSummary_5min_meta.csv | LiP-MS data; 5 min after refolding initiated, sample 2 (reformatted from matched .xlsx) | Fried Lab, unpublished |
| 12 | 20220505_YeastRefolding_ProteinSummary_5min_duplicate_meta.xlsx | LiP-MS data; 5 min after refolding initiated, sample 3 | Fried Lab, unpublished |
| 12 | 20220505_YeastRefolding_ProteinSummary_5min_duplicate_meta.csv |  LiP-MS data; 5 min after refolding initiated, sample 3 (reformatted from matched .xlsx) | Fried Lab, unpublished |
| 12 | 20210721_YeastRefold_ProteinSummary_2hr_meta.xlsx | LiP-MS data; 120 min after refolding initiated, sample 1 | Fried Lab, unpublished |
| 12 | 20210721_YeastRefold_ProteinSummary_2hr_meta.csv | LiP-MS data; 120 min after refolding initiated, sample 1 (reformatted from matched .xlsx) | Fried Lab, unpublished |
| 12 | 20220505_YeastRefolding_ProteinSummary_2hr_meta.xlsx | LiP-MS data; 120 min after refolding initiated, sample 2 | Fried Lab, unpublished |
| 12 | 20220505_YeastRefolding_ProteinSummary_2hr_meta.csv | LiP-MS data; 120 min after refolding initiated, sample 2 (reformatted from matched .xlsx) | Fried Lab, unpublished |
| 12 | 20220314_HeatShock_ProteinSummary_meta.xlsx | LiP-MS data, heatshocked versus normal cells | Fried Lab, unpublished |
| 12 | 20220314_HeatShock_ProteinSummary_meta.csv | LiP-MS data, heatshocked versus normal cells (reformatted from matched .xlsx) | Fried Lab, unpublished |
| 12 | 20220314_Recovery_ProteinSummary_meta.xlsx | LiP-MS data, normal cells versus cells recovering from heatshock | Fried Lab, unpublished | 
| 12 | 20220314_Recovery_ProteinSummary_meta.csv | LiP-MS data, normal cells versus cells recovering from heatshock (reformatted from matched .xlsx) | Fried Lab, unpublished |
| 13 | Yeast_AF_combined_20250530.csv | Entanglement status of proteins from AlphaFold2 structures | O'Brien Lab, unpublished |
| 13 | Yeast_EXP_combined_20250530.csv | Entanglement status of proteins from experimental structures | O'Brien Lab, unpublished |
| 14 | 1-s2.0-S2211124717312160-mmc2.tsv | Chaperones and cochaperones in S. cerevisiae| 10.1016/j.celrep.2017.08.074, Table S1 |
| 15 | 559292.tsv | Complex Portal S. cerevisiae protein complexes | Downloaded from https://ftp.ebi.ac.uk/pub/databases/intact/complex/current/complextab/559292.tsv |
| 17 | inviable_annotations.txt | Essential yeast proteins from Saccharomyces Genome Database | Downloaded from yeastgenome.org |
| 17 | inviable_annotations-mod.txt | Essential yeast protein from Saccharomyces Genome Database, reformatted | Downloaded from yeastgenome.org |
| 18 | Y2H_union.txt | Yeast two-hybrid protein-protein interaction data | Downloaded from https://interactome.dfci.harvard.edu/S_cerevisiae/download/Y2H_union.txt |
| 18 | Y2H_union-clean.txt | Yeast two-hybrid data with protein names cleaned | See the program `18-Y2H-data/python-scripts/clean-Y2H-union.py` |
| 18 | Y2H_union-clean_edges.csv | List of edges (interactions) between nodes | Extracted using Cytoscape |
| 18 | Y2H_union-clean_nodes.csv | List of nodes (proteins) | Extracted using Cytoscape |
| 19 | meltome-atlas-ma_0010.txt | Meltome Atlas data for S. cerevisiae proteome | 10.1038/s41592-020-0801-4, extracted from Table S2 |

