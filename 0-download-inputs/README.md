### `0-download-inputs`: Download input files used by downstream pipeline steps

Run the command:

```bash
python check-node-count.py --file1 data-files/ppi_in_vivo_edges.csv --file2 data-files/ppi_in_vitro_edges.csv --merged_edges data-files/merged-edges.csv --merged_nodes data-files/merged-nodes.csv
```

To generate the merged node and edge files



Corresponding configuration file section: `DOWNLOAD INPUT FILES`

Anticipated execution time: 60-120 min

This section of the overall pipeline handles the download of required files to `0-download-inputs/data-files/`

* To run this step in isolation, navigate to the main project directory (i.e., up one level from `0-download-inputs`) and run the command:

```bash
snakemake -c all --use-conda --conda-frontend conda --snakefile 0-download-inputs/Snakefile --configfile config-files/s288c.config
```

* The runtime of this pipeline depends on connection speeds

#### Requirements for running this pipeline

This pipeline requires that have a working `gocmd` executable on your machine. You will need:

1. a CyVerse account (sign up [here](https://user.cyverse.org/register))
2. the `gocmd` executable (download it [here](https://learning.cyverse.org/ds/gocommands/installation/))
3. to `upgrade` your executable with `sudo gocmd upgrade`
4. initialize 'gocmd' with `gocmd init` (input your credentials as below)

`gocmd init` credentials:

Upon running this command, you will be prompted to input five pieces of information in series on the command line:
  * (1) `iRODS Host [data.cyverse.org]`: hit enter to accept the default of `data.cyverse.org`
  * (2) `iRODS Port [1247]`: hit enter to accept the default of `1247`
  * (3) `iRODS Zone [iplant]`: hit enter to accept the default of `iplant`
  * (4) `iRODS Username`: type in your CyVerse username and then hit enter
  * (5) `iRODS Password`: type in your CyVerse password and then hit enter. Note well: you will not see your password or asterisks representing the characters you have entered, but your keystrokes are being recorded. 

You are now ready to run this phase of the overall pipeline - be sure to insert the path to your `gocmd` executable in the `config` file.

#### Summary of data files

No input data files are distributed on GitHub. All data for *S. cerevisiae* are archived on CyVerse Data Store at the path 
```text
/iplant/home/shared/NCEMS/working-groups/energetic-origins/arabidopsis-thaliana/0-download-inputs/data-files`
```

This path is encoded in `../config-files/athaliana.config`

#### Current data files downloaded

* `arabidopsis_candidate_pairs_ranked.tsv`: the post-processed interactions from DOI: 10.1016/j.cell.2020.02.049 (supplementary files S4 and https://doi.org/10.5281/zenodo.3666940)

* currently the pipeline gives an error message because the file download counts are incorrect

### Legacy README.md from s288c

The following 31 files will be downloaded:

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
|20|13| Yeast_AF_combined_20250530.csv | Entanglement status of proteins from AlphaFold2 structures | O'Brien Lab, unpublished |
|21|13| Yeast_EXP_combined_20250530.csv | Entanglement status of proteins from experimental structures | O'Brien Lab, unpublished |
|22|14| 1-s2.0-S2211124717312160-<br>mmc2.tsv | Chaperones and cochaperones in S. cerevisiae | 10.1016/j.celrep.2017.08.074, Table S1 |
|23|15| 559292.tsv | Complex Portal S. cerevisiae protein complexes | Downloaded from [Complex Portal](https://ftp.ebi.ac.uk/pub/databases/intact/complex/current/complextab/559292.tsv) |
|24|15| YEAST_559292_idmapping.dat | Yeast ID mapping file from UniProt | Downloaded from [UniProt](https://ftp.uniprot.org/pub/databases/uniprot/knowledgebase/idmapping/by_organism/YEAST_559292_idmapping.dat.gz) |
|25|17| inviable_annotations.txt | Essential yeast proteins from Saccharomyces Genome Database | Downloaded from yeastgenome.org |
|26|17| inviable_annotations-mod.txt | Essential yeast protein from Saccharomyces Genome Database, reformatted | Downloaded from yeastgenome.org |
|27|18| Y2H_union.txt | Yeast two-hybrid protein-protein interaction data | [Downloaded from CCSB Interactome Database](https://interactome.dfci.harvard.edu/S_cerevisiae/download/Y2H_union.txt) |
|28|18| Y2H_union-clean.txt | Yeast two-hybrid data with protein names cleaned | See the program `18-Y2H-data/python-scripts/clean-Y2H-union.py` |
|29|18| Y2H_union-clean_edges.csv | List of edges (interactions) between nodes | Extracted using Cytoscape |
|30|18| Y2H_union-clean_nodes.csv | List of nodes (proteins) | Extracted using Cytoscape |
|31|19| meltome-atlas-ma_0010.txt | Meltome Atlas data for S. cerevisiae proteome | 10.1038/s41592-020-0801-4, extracted from Table S2 |

* After running the pipeline (downloading and unpacking), you should have the following file counts in `0-download-inputs/data-files`:

| File extension | Number of files |
|:--------------:|----------------:|
| fasta          | 6,040           |
| pdb            | 6,039           |
| csv            |    11           |
| txt            |     5           |
| tsv            |     4           |
| xlsx           |     3           |
| dat            |     2           |
| obo            |     1           |
| pt             |     1           |
| xml            |     1           |
| cys            |     1           |
| TOTAL FILES    |12,108           |

These file counts are verified by the pipeline at the end of the run. 

* Once pipeline step `2-sequence-parsing` is executed, an additional 129 PDB files will be added for truncated AlphaFold2 PDBs with signal sequences removed. This brings the number of PDB files to 6,168 and the total number of files to 12,237.
