# Figure 4 — Predictive Features of Network Hubs (Random Forest)

This folder contains the Jupyter notebooks used to generate **Figure 4** of the manuscript. These notebooks implement **Random Forest classification models** to identify features associated with hub status in the yeast interaction networks.

## Overview

We trained **four supervised random forest classification models**, each corresponding to a different network and centrality definition:

1. **AE-MS Degree Centrality**
2. **AE-MS Betweenness Centrality**
3. **Y2H Degree Centrality**
4. **Y2H Betweenness Centrality**

Instead of a single train/test split, each model was evaluated using an **ensemble-style repeated training framework** of **50 Random Forest runs**. In this approach, the data were repeatedly split using stratified shuffle splitting, and a new Random Forest model was trained on each split. Final performance was summarized across the 50 runs using the **mean ± standard deviation** of classification metrics.

This repeated-model framework was used to provide a more robust estimate of model performance and feature importance than a single split alone.

Each notebook outputs:

- an **Excel file** containing the ROC curve coordinates (`False_Positive_Rate`, `True_Positive_Rate`)
- a beeswarm SHAP plot

These four ROC Excel files are then combined to generate **Figure 4A** of the manuscript.













