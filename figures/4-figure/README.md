# Figure 4 — Predictive Features of Network Hubs (Random Forest)

This folder contains Jupyter notebooks for four supervised classification models using Random Forests.  
Each model has four variants (v1–v4, v5–v8, v9–v12, v13–v16) corresponding to different training strategies.  
For each task, the **best-performing model** (marked **FINAL**) was selected for Figure 4.

- **Figure 4A**: ROC curves for the four final models (MS-degree, MS-betweenness, Y2H-degree, Y2H-betweenness), generated from their `roc_model*.csv` files.  
- **Figure 4B**: SHAP summary for the **MS-degree** final model (`SHAP_random_forest3.svg`).  
- **Figure 4C**: SHAP summary for the **Y2H-degree** final model (`SHAP_random_forest12.svg`). 

---

## Tasks & Notebooks

### 1) Predict **MS-degree centrality**
**Variants (v1–v4):**
| Variant | Training strategy              | Notebook filename                                                |
|:-------:|--------------------------------|------------------------------------------------------------------|
| v1      | No Hyperparameter Tuning       | `v-4.1-RF-Classifier-MS-degree-percentile.ipynb`                 |
| v2      | RandomizedSearchCV             | `v-4.2-RF-Classifier-MS-degree-percentile.ipynb`                 |
| v3      | **GridSearchCV (FINAL)**       | `v-4.3-RF-Classifier-MS-degree-percentile-FINAL.ipynb`           |
| v4      | SMOTE + GridSearchCV           | `v-4.4-RF-Classifier-MS-degree-percentile.ipynb`                 |

**Final outputs:**
- ROC: `roc_model3.csv` (contains FPR, TPR)
- SHAP: `SHAP_random_forest1.svg, SHAP_random_forest2.svg, SHAP_random_forest3.svg, SHAP_random_forest4.svg`

---

### 2) Predict **MS-betweenness centrality**

**Variants (v5–v8):**
| Variant | Training strategy              | Notebook filename                                                |
|:-------:|--------------------------------|------------------------------------------------------------------|
| v5      | No Hyperparameter Tuning       | `v-4.5-RF-Classifier-MS-degree-percentile.ipynb`                 |
| v6      | **RandomizedSearchCV (FINAL)** | `v-4.6-RF-Classifier-MS-degree-percentile-FINAL.ipynb`           |
| v7      | GridSearchCV                   | `v-4.7-RF-Classifier-MS-degree-percentile.ipynb`                 |
| v8      | SMOTE + GridSearchCV           | `v-4.8-RF-Classifier-MS-degree-percentile.ipynb`                 |

**Final outputs:**
- ROC: `roc_model6.csv` (contains FPR, TPR)
- SHAP: `SHAP_random_forest5.svg, SHAP_random_forest6.svg, SHAP_random_forest7.svg, SHAP_random_forest8.svg`

---

### 3) Predict **Y2H-degree centrality**
**Variants (v9–v12):**
| Variant | Training strategy              | Notebook filename                                                |
|:-------:|--------------------------------|------------------------------------------------------------------|
| v9      | No Hyperparameter Tuning       | `v-4.9-RF-Classifier-Y2H-degree-percentile.ipynb`                |
| v10     | RandomizedSearchCV             | `v-4.10-RF-Classifier-Y2H-degree-percentile.ipynb`               |
| v11     | GridSearchCV                   | `v-4.11-RF-Classifier-Y2H-degree-percentile.ipynb`               |
| v12     | **SMOTE + GridSearchCV (FINAL)** | `v-4.12-RF-Classifier-Y2H-degree-percentile-FINAL.ipynb`       |

**Final outputs:**
- ROC: `roc_model12.csv` (contains FPR, TPR)
- SHAP: `SHAP_random_forest9.svg, SHAP_random_forest10.svg, SHAP_random_forest11.svg, SHAP_random_forest12.svg`

---

### 4) Predict **Y2H-betweenness centrality**
> *Note:* Original filename spelling retained for v14: “betwenness”.

**Variants (v13–v16):**
| Variant | Training strategy              | Notebook filename                                                |
|:-------:|--------------------------------|------------------------------------------------------------------|
| v13     | No Hyperparameter Tuning       | `v-4.13-RF-Classifier-Y2H-betweenness-percentile.ipynb`          |
| v14     | RandomizedSearchCV             | `v-4.14-RF-Classifier-Y2H-betwenness-percentile.ipynb`           |
| v15     | **GridSearchCV (FINAL)**       | `v-4.15-RF-Classifier-Y2H-betweenness-percentile-FINAL.ipynb`    |
| v16     | SMOTE + GridSearchCV           | `v-4.16-RF-Classifier-Y2H-betweenness-percentile.ipynb`          |

**Final outputs:**
- ROC: `roc_model15.csv` (contains FPR, TPR)
- SHAP: SHAP: `SHAP_random_forest13.svg, SHAP_random_forest14.svg, SHAP_random_forest15.svg, SHAP_random_forest16.svg`

---

## Outputs

- **ROC CSVs (used for Figure 4A):**
  - `roc_model3.csv` — MS-degree (Task 1 FINAL)
  - `roc_model6.csv` — MS-betweenness (Task 2 FINAL)
  - `roc_model12.csv` — Y2H-degree (Task 3 FINAL)
  - `roc_model15.csv` — Y2H-betweenness (Task 4 FINAL)

- **SHAP summaries (one per variant):**
  - `SHAP_random_forest1.svg` … `SHAP_random_forest16.svg`  
    - Mapping:  
      - 1–4 → Task 1 (v1–v4), **3** is FINAL  
      - 5–8 → Task 2 (v5–v8), **6** is FINAL  
      - 9–12 → Task 3 (v9–v12), **12** is FINAL  
      - 13–16 → Task 4 (v13–v16), **15** is FINAL

---


