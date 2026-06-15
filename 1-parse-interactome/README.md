# 1-parse-interactome

This pipeline step curates, harmonizes, and reformats Arabidopsis thaliana protein-protein interaction edge lists for downstream network annotation and analysis.

## Purpose

The goal of this step is to generate clean in vivo, in vitro, union, and intersection protein interaction networks using a consistent identifier namespace. Input interaction data may contain a mixture of AGI locus identifiers, UniProt accessions, and KEGG-style identifiers. Before network construction, all edge endpoints are mapped to Arabidopsis AGI locus IDs whenever possible.

Using a single identifier system at this stage is critical because downstream steps rely on AGI IDs for sequence assignment, isoform selection, protein annotation, and integration with other Arabidopsis datasets.

## Inputs

The pipeline expects the following input files, defined in the project config file:

- `BIOGRID-ORGANISM-Arabidopsis_thaliana_Columbia-5.0.258.tab3.xlsx`
  - BIOGRID interaction data for Arabidopsis thaliana.
- Manually annotated Yilmaz/Table S1 interaction file
  - Expected to contain method/evidence annotations used to classify interactions.
- `ARATH_3702_idmapping.dat`
  - UniProt Arabidopsis identifier mapping file.
  - Used to map UniProt and KEGG-style identifiers to AGI locus IDs.

## Main processing steps

### 1. Curate raw interaction evidence

The `process_data` rule runs:

```bash
curate-athaliana-ppi-evidence.py
```

This script reads BIOGRID and the manually annotated Table S1 interaction file, classifies interactions as in vivo, in vitro, or unknown, and writes raw curated edge files:

```text
ppi-in-vivo-edges.raw.csv
ppi-in-vitro-edges.raw.csv
ppi-unknown-methods.csv
```

### 2. Map identifiers to AGI locus IDs

The map_identifiers rule runs:

```bash
map-edge-identifiers-to-agi.py
```

This script maps all edge endpoints to AGI locus IDs using the UniProt idmapping.dat file.

Supported input identifier types include:

AGI locus IDs, e.g. AT1G01010
AGI isoform-like IDs, e.g. AT1G01010.1
UniProt accessions, e.g. Q9FEN9
UniProt isoform accessions, e.g. Q8VZJ1-1
KEGG-style identifiers, e.g. ath:ArthCp058 or ArthCp058

Mappings are conservative. Identifiers are retained only if they are already AGI IDs or can be mapped unambiguously to a single AGI ID. Edges with unmapped, missing, or ambiguous endpoints are dropped.

This rule writes:

```text
ppi-in-vivo-edges.csv
ppi-in-vitro-edges.csv
identifier-mapping-in-vivo-report.csv
identifier-mapping-in-vitro-report.csv
```

The mapping reports should be inspected to evaluate dropped or ambiguous identifiers.

3. Build union and intersection networks

The merge_and_intersect rule runs:

```bash
merge-and-intersect.py
```

This script creates reformatted individual edge lists, individual node lists, union network files, and node-intersection-filtered edge lists.

Outputs include:

```text
union-nodes.csv
union-edges.csv
in-vivo-nodes.csv
in-vivo-edges.csv
in-vitro-nodes.csv
in-vitro-edges.csv
intersection-nodes.csv
intersection-in-vivo-edges.csv
intersection-in-vitro-edges.csv
```

The intersection node set is defined as the set of AGI identifiers present in both the in vivo and in vitro networks. The intersection edge files are not edge-set intersections; instead, they contain the original in vivo or in vitro edges filtered to edges whose endpoints are both present in the shared node set.

Outputs

The primary downstream outputs are:

| File | Description |
| :----------: | :----------: |
| union-nodes.csv | Node list for the union of in vivo and in vitro networks |
| union-edges.csv | Edge list for the union network |
| in-vivo-nodes.csv | Node list for the AGI-normalized in vivo network |
| in-vivo-edges.csv | Reformatted AGI-normalized in vivo edge list |
| in-vitro-nodes.csv | Node list for the AGI-normalized in vitro network |
| in-vitro-edges.csv | Reformatted AGI-normalized in vitro edge list |
| intersection-nodes.csv | Nodes present in both in vivo and in vitro networks |
| intersection-in-vivo-edges.csv | In vivo edges possible within the intersection node set |
| intersection-in-vitro-edges.csv | In vitro edges possible within the intersection node set |
| identifier-mapping-in-vivo-report.csv | Mapping report for in vivo edge identifiers |
| identifier-mapping-in-vitro-report.csv | Mapping report for in vitro edge identifiers |

### Running the workflow

From the repository root:

```bash
snakemake -c all \
  --use-conda \
  --conda-frontend conda \
  --configfile config-files/union-athaliana.config \
  --snakefile 1-parse-interactome/Snakefile
```

### Quality checks

After the workflow completes, confirm that the final node list contains only AGI-like identifiers:

awk -F, 'NR > 1 && $1 !~ /^AT[1-5CM]G[0-9]{5}$/ {print $1}' \
  1-parse-interactome/processed-data/union-nodes.csv | head

The expected result is no output.

Count any remaining non-AGI nodes:

awk -F, 'NR > 1 && $1 !~ /^AT[1-5CM]G[0-9]{5}$/ {n++} END {print n+0}' \
  1-parse-interactome/processed-data/union-nodes.csv

The expected result is:

0

Inspect unresolved identifier mappings:

grep ",unmapped" 1-parse-interactome/processed-data/identifier-mapping-*-report.csv | head
grep ",ambiguous" 1-parse-interactome/processed-data/identifier-mapping-*-report.csv | head
Notes

Generated files in processed-data/ should generally not be committed to Git. They can be regenerated by Snakemake from the raw input files and scripts.

The identifier mapping step intentionally favors specificity over recall. Some edges may be dropped if one or both endpoints cannot be confidently mapped to a unique AGI locus ID. This is preferable to carrying mixed identifier types into downstream analyses.
