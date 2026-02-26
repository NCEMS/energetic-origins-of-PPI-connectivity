This directory contains Jupyter notebooks used to generate all main text and supplementary figures for the manuscript.
Each subfolder corresponds to figures and contains the notebook(s) that render the final .svg files.

  **Directory Overview**

1. **1-figure/** → Figure 1. Construction and Characterization of Network Models for the S. cerevisiae Interactome.
2. **2-figure/** → Figure 2. Cellular Abundance and Turnover of Network Hubs and Bottlenecks.
3. **3-figure/** → Figure 3. Conformational Stability of Hub Proteins.
4. **4-figure/** → Figure 4. Predictive Features of Network Hubs.
   - Model used: Random Forest Classiffier
   - Feature sets: 12 feature sets defined in the manuscript
   - Evaluation: AUROC, AUPRC, accuracy, F1, and cross-validation
   - Reproducibility: Random seeds and cross-validation strategy are set within the notebooks
   - Contains multiple notebooks that train several ML models and generate all Fig. 4 panels.
   - Refer to the README.md inside this folder for detailed structure.
5. **5-figure/** → Figure 5. Context-Dependent Properties of Hub Proteins.
6. **6-figure/** → Figure 6. Chaperone Dependence of Hub Proteins.
7. **supplementary/** → Supplementary Figures S1–S7.

All generated figures appear in their respective folders, saved as `.svg` files with filenames matching the manuscript figure numbering/panel order.


 **Data Inputs**

- All notebooks contain explicit references to the data files they load.
- Data paths and preprocessing steps are documented within each notebook.
- There are no external preprocessing scripts required to run any figure notebook.

 **Running the Notebooks**

- No required run order.
- No cross-notebook dependencies.
- All notebooks can be executed independently.

 **Outputs**

- Each notebook produces all figure panels within the same folder where the notebook resides.
- All final figure files:
  - are saved as SVG  
  - follow the manuscript figure numbering  
  - appear in the correct panel order  
