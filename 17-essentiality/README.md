### `17-essentiality`: Integrate information about essential proteins

* This pipeline merges proteins annotated within the Saccharomyces Genome Database as essential ("inviable" within SGD)
* Data were downloaded from SGD by navigating through the following menu selections: select "Function" > select "Phenotype" > select "Browse All Phenotypes" > select "Inviable" > download all entries

Note: 

The file downloaded from SGD (`../0-download-inputs/data-files/inviable_annotations.txt`) has several malformed lines. 

The following manual modifications were made to this file:

Lines 1558 & 1559 were merged by deleting the newline after line 1558

Lines 1560 & 1561 were merged by deleting the newline after line 1560

The result is `../0-download-inputs/data-files/inviable_annotations-mod.txt`, which is read into the pipeline in this step.
