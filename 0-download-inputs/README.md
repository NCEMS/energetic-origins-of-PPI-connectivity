### `0-download-inputs`: Download input files used by downstream pipeline steps

* This section of the overall pipeline handles the download of required files to `0-download-inputs/data-files/`

* To run this step in isolation, navigate to the main project directory (i.e., up one level from `0-download-inputs`) and run the command:

`snakemake --snakefile 0-download-inputs/Snakefile --configfile config-files/s288c.config --c all --use-conda --conda-frontend conda`

* This pipeline will require 60-90 min depending on connection speeds and the write speed of your file system. 

*N.B.*: Once all files are downloaded and unpacked, `0-download-inputs/data-files` will contain ~112 GB of data. 

#### Summary of data files

The following 43 files are distributed with the GitHub repository. In most cases, this is because the file required pre-processing not easily accomplished programmatically (e.g., converting from .xlsx to .csv). In such cases, the original file is also included in the repository (see, for example, `1-s2.0-S2211124714009346-mmc2.xlsx` and `1-s2.0-S2211124714009346-mmc2.csv`)

|Index |Used in step| Filename | Description | Source |
|:----:|:----------:|:---------|:------------|:-------|
|1|1| The_Yeast_Interactome.cys | The Yeast Interactome Cytoscape session file | Downloaded from [The Yeast Interactome](https://www.yeast-interactome.org/) |
|2|1| The_Yeast_Interactome_nodes.csv | List of nodes (proteins) | Extracted from .cys file with Cytoscape |
|3|1| The_Yeast_Interactome_edges.csv | List of edges (interactions) between nodes | Extracted from .cys file with Cytoscape |
|4|2| orf_trans.fasta | S. cerevisiae protein sequences from Saccharomyces Genome Database| [SGD Archive](http://sgd-archive.yeastgenome.org/sequence/S288C_reference/orf_protein/orf_trans.fasta.gz)|
|5|3| uniprot_sprot.xml | UniProt database in .xml format| [UniProt](https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.xml.gz) |
|6|3| go-basic.obo| Open Biomedical Ontologies data file | [OboFoundry](http://purl.obolibrary.org/obo/go.obo) |
|7|3| Yeast_GSB_phospho_all_prots_0125.csv | ProteomeeXchange S. cerevisiae post-translational modifications| [ProteomeeXchange](https://www.proteomexchange.org/ptmexchange/Files/Yeast_GSB_phospho_all_prots_0125.csv.gz) |
|8|4| DisProt-release_2024_12-<br>with_ambiguous_evidences.tsv | DisProt database version 2024_12 | Downloaded from [DisProt](https://disprot.org/download) |
|9|4| DisProt-release_2024_12-<br>with_ambiguous_evidences-<br>cleaned.tsv | DisProt database version 2024_12 | Cleaned version of DisProt with whitespace cleaned for reading; generated with `bash-scripts/clean-disprot.sh` outside of main pipeline |
|10|5| UP000002311_559292_YEAST_v4.tar | S. cerevisiae AlphaFold2 structure predictions| [EBI](https://ftp.ebi.ac.uk/pub/databases/alphafold/latest/UP000002311_559292_YEAST_v4.tar)|
|11|5| esm_if1_gvp4_t16_142M_UR50.pt | ESM-IF model for Cagiada et al. 2025 stability predictions | [Electronic Research Data Archive](https://sid.erda.dk/share_redirect/eIZVVNEd8B) |
|12|6| protein2ipr.dat| InterPro domain annotations | [InterPro](https://ftp.ebi.ac.uk/pub/databases/interpro/releases/latest/protein2ipr.dat.gz) |
|13|8| 1-s2.0-S2211124714009346-<br>mmc2.xlsx | Protein half-lives in S. cerevisiae | 10.1016/j.celrep.2014.10.065, Table S1 |
|14|8| 1-s2.0-S2211124714009346-<br>mmc2.csv | Protein half-lives in S. cerevisiae | 10.1016/j.celrep.2014.10.065, Table S1 (reformatted for analysis) |
|15|8| 1-s2.0-S2405471217303411-<br>mmc2.xlsx | Protein half-lives in S. cerevisiae | 10.1016/j.cels.2017.08.008, Data S1 |
|16|8| 1-s2.0-S2405471217303411-<br>mmc2.csv | Protein half-lives in S. cerevisiae | 10.1016/j.cels.2017.08.008, Data S1 (reformatted for analysis) |
|17|9| 1-s2.0-S240547121730546X-<br>mmc5.xlsx | Protein expression in S. cerevisiae | 10.1016/j.cels.2017.12.004, Table S4 |
|18|9| 1-s2.0-S240547121730546X-<br>mmc5.csv | Protein expression in S. cerevisiae | 10.1016/j.cels.2017.12.004, Table S4 (reformatted for analysis)|
|19|10| skr_weinberg_genesTE.csv | Translation efficiency information from scikit-ribo for S. cerevisiae | Downloaded from [GitHub](https://github.com/schatzlab/scikit-ribo_manuscript/blob/master/Data/skr_weinberg_genesTE.csv) |
|20|12| 20210721_YeastRefold_ProteinSu<br>mmary_1min_meta.xlsx | LiP-MS data; 1 min after refolding initiated, sample 1 | Fried Lab, unpublished |
|21|12| 20210721_YeastRefold_ProteinSu<br>mmary_1min_meta.csv | LiP-MS data; 1 min after refolding initiated, sample 1 (reformatted from matched .xlsx) | Fried Lab, unpublished |
|22|12| 20220505_YeastRefolding_Protei<br>nSummary_1min_meta.xlsx | LiP-MS data; 1 min after refolding initiated, sample 2 | Fried Lab, unpublished |
|23|12| 20220505_YeastRefolding_Protei<br>nSummary_1min_meta.csv | LiP-MS data; 1 min after refolding initiated, sample 2 (reformatted from matched .xlsx) | Fried Lab, unpublished |
|24|12| 20220505_YeastRefolding_Protei<br>nSummary_1min_duplicate_meta.xlsx | LiP-MS data; 1 min after refolding initiated, sample 3 | Fried Lab, unpublished |
|25|12| 20220505_YeastRefolding_Protei<br>nSummary_1min_duplicate_meta.csv | LiP-MS data; 1 min after refolding initiated, sample 3 (reformatted from matched .xlsx) | Fried Lab, unpublished |
|26|12| 20210721_YeastRefold_ProteinSu<br>mmary_5min_meta.xlsx | LiP-MS data; 5 min after refolding initiated, sample 1 | Fried Lab, unpublished |
|27|12| 20210721_YeastRefold_ProteinSu<br>mmary_5min_meta.csv | LiP-MS data; 5 min after refolding initiated, sample 1 (reformatted from matched .xlsx) | Fried Lab, unpublished |
|28|12| 20220505_YeastRefolding_Protei<br>nSummary_5min_meta.xlsx | LiP-MS data; 5 min after refolding initiated, sample 2 | Fried Lab, unpublished |
|29|12| 20220505_YeastRefolding_Protei<br>nSummary_5min_meta.csv | LiP-MS data; 5 min after refolding initiated, sample 2 (reformatted from matched .xlsx) | Fried Lab, unpublished |
|30|12| 20220505_YeastRefolding_Protei<br>nSummary_5min_duplicate_meta.xlsx | LiP-MS data; 5 min after refolding initiated, sample 3 | Fried Lab, unpublished |
|31|12| 20220505_YeastRefolding_Protei<br>nSummary_5min_duplicate_meta.csv | LiP-MS data; 5 min after refolding initiated, sample 3 (reformatted from matched .xlsx) | Fried Lab, unpublished |
|32|12| 20210721_YeastRefold_ProteinSu<br>mmary_2hr_meta.xlsx | LiP-MS data; 120 min after refolding initiated, sample 1 | Fried Lab, unpublished |
|33|12| 20210721_YeastRefold_ProteinSu<br>mmary_2hr_meta.csv | LiP-MS data; 120 min after refolding initiated, sample 1 (reformatted from matched .xlsx) | Fried Lab, unpublished |
|34|12| 20220505_YeastRefolding_Protei<br>nSummary_2hr_meta.xlsx | LiP-MS data; 120 min after refolding initiated, sample 2 | Fried Lab, unpublished |
|35|12| 20220505_YeastRefolding_Protei<br>nSummary_2hr_meta.csv | LiP-MS data; 120 min after refolding initiated, sample 2 (reformatted from matched .xlsx) | Fried Lab, unpublished |
|36|12| 20220314_HeatShock_ProteinSumm<br>ary_meta.xlsx | LiP-MS data, heatshocked versus normal cells | Fried Lab, unpublished |
|37|12| 20220314_HeatShock_ProteinSumm<br>ary_meta.csv | LiP-MS data, heatshocked versus normal cells (reformatted from matched .xlsx) | Fried Lab, unpublished |
|38|12| 20220314_Recovery_ProteinSummar<br>y_meta.xlsx | LiP-MS data, normal cells versus cells recovering from heatshock | Fried Lab, unpublished |
|39|12| 20220314_Recovery_ProteinSummar<br>y_meta.csv | LiP-MS data, normal cells versus cells recovering from heatshock (reformatted from matched .xlsx) | Fried Lab, unpublished |
|40|13| Yeast_AF_combined_20250530.csv | Entanglement status of proteins from AlphaFold2 structures | O'Brien Lab, unpublished |
|41|13| Yeast_EXP_combined_20250530.csv | Entanglement status of proteins from experimental structures | O'Brien Lab, unpublished |
|42|14| 1-s2.0-S2211124717312160-<br>mmc2.tsv | Chaperones and cochaperones in S. cerevisiae | 10.1016/j.celrep.2017.08.074, Table S1 |
|43|15| 559292.tsv | Complex Portal S. cerevisiae protein complexes | Downloaded from [Complex Portal](https://ftp.ebi.ac.uk/pub/databases/intact/complex/current/complextab/559292.tsv) |
|44|17| inviable_annotations.txt | Essential yeast proteins from Saccharomyces Genome Database | Downloaded from yeastgenome.org |
|45|17| inviable_annotations-mod.txt | Essential yeast protein from Saccharomyces Genome Database, reformatted | Downloaded from yeastgenome.org |
|46|18| Y2H_union.txt | Yeast two-hybrid protein-protein interaction data | [Downloaded from CCSB Interactome Database](https://interactome.dfci.harvard.edu/S_cerevisiae/download/Y2H_union.txt) |
|47|18| Y2H_union-clean.txt | Yeast two-hybrid data with protein names cleaned | See the program `18-Y2H-data/python-scripts/clean-Y2H-union.py` |
|48|18| Y2H_union-clean_edges.csv | List of edges (interactions) between nodes | Extracted using Cytoscape |
|49|18| Y2H_union-clean_nodes.csv | List of nodes (proteins) | Extracted using Cytoscape |
|50|19| meltome-atlas-ma_0010.txt | Meltome Atlas data for S. cerevisiae proteome | 10.1038/s41592-020-0801-4, extracted from Table S2 |

* After running this pipeline step, you should have the following file counts in `0-download-inputs/data-files`:

| File extension | Number of files |
|:--------------:|----------------:|
| pdb            | 6,039           |
| fasta          | 6,040           |
| csv            |    21           |
| xlsx           |    13           |
| txt            |     5           |
| tsv            |     4           |
| dat            |     2           |
| obo            |     1           |
| pt             |     1           |
| xml            |     1           |
| cys            |     1           |
| TOTAL FILES    |12,128           |

* Once pipeline step `2-sequence-parsing` is executed, an additional 129 PDB files will be added corresponding to truncated AlphaFold2 PDBs with signal sequences removed. This brings the number of PDB files to 6,168 and the total number of files to 12,257.
