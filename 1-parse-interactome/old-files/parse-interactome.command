python python-scripts/merge-and-intersect.py \
--file1 data-files/ppi_in_vitro_edges.csv \
--file2 data-files/ppi_in_vivo_edges.csv \
--union_edges data-files/merged-edges.csv \
--union_nodes data-files/merged-nodes.csv \
--file1_nodes data-files/in-vitro-nodes.csv \
--file2_nodes data-files/in-vivo-nodes.csv \
--file1_intersection_edges data-files/intersection-edges-in_vitro.csv \
--file2_intersection_edges data-files/intersection-edges-in_vivo.csv \
--intersection_nodes data-files/intersection-nodes.csv
