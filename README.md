# 🧬 scIBD Colon – Cluster-wise Gene Expression Summary & Disease-specific Significance

This repository provides a reproducible pipeline for analyzing **single-cell colon data (scIBD)** to identify genes whose expression levels differ across diseases and within cell clusters.  

The workflow:
1. **Cell counts per disease**  
   Produces: `outputs/cell_type_count_per_disease.csv`

2. **Per-disease gene expression statistics per cell type**  
   For every cell type, computes **min, max, mean, variance** of every gene for every disease.  
   Produces: `outputs/stats_per_celltype/{celltype}_stats.csv`

3. **Kruskal–Wallis test per gene per cell type**  
   Non-parametric multi-group test across diseases.  
   Produces: `outputs/kw_per_celltype/{celltype}_KW.csv`

4. **Effect size per gene per cell type**  
   Effect size defined as:  
   \[
   \text{effect\_size} = \max(\text{mean across diseases}) - \min(\text{mean across diseases})
   \]
   and categorized into `low` (<0.25), `medium` (0.25–0.5), `high` (>0.5).  
   Produces: `outputs/effectsize_per_celltype/{celltype}_effectsize.csv`

5. **Pairwise post-hoc testing**  
   For every gene, every cell type, and every disease pair, runs Mann–Whitney U (two-sided) and adjusts p-values with Benjamini–Hochberg.  
   Produces: `outputs/posthoc_per_celltype/{celltype}_pairwise_posthoc.csv`

6. **UC_inflamed-specific separation**  
   Selects genes *per cell type* for which:
   - Kruskal–Wallis p < 0.05 / number of genes
   - Kruskal–Wallis FDR < 0.05
   - Effect size > 0.25
   - All pairwise tests involving `UC_inflamed` pass
     - Mann–Whitney p < 0.05 / number of genes
     - Mann–Whitney FDR < 0.05  
   Produces: `outputs/uc_inflamed_separated/{celltype}_UCinflamed.csv`

---

## ❓ Research Question

> *Which genes exhibit disease-associated expression differences in single-cell colon data, and how consistent are these effects across cell clusters?*

---

## 📥 Input

- **AnnData file:**  
  `./scIBD_Colon.h5ad`

- **Required columns in** `adata.obs`  
  - `disease`
  - `major_cluster`

- **Gene names:** stored in `adata.var_names`  
- **Expression matrix:** `adata.X` (can be sparse or dense)

---

## ⚙️ Requirements

Python ≥ 3.9  
Install dependencies:
```bash
pip install scanpy pandas numpy scipy
````

Make sure to include this import in your script:

```python
from scipy import sparse as sp
```

---

## 🚀 How to Run

```bash
python main.py
```

This script:

1. Loads the `.h5ad` file using **Scanpy**.
2. Computes per-cluster gene statistics.
3. Runs Kruskal–Wallis tests.
4. calculates effect size per gene per cell type.
6. Runs pairwise post-hoc testin.
7. Evaluates the results for UC_inflamed disease.

---

## 📂 Output Structure

```
.
├── data/
│   └── input.h5ad
├── outputs/
│   ├── cell_type_count_per_disease.csv
│   ├── stats_per_celltype/
│   ├── kw_per_celltype/
│   ├── effectsize_per_celltype/
│   ├── posthoc_per_celltype/
│   └── uc_inflamed_separated/
└── main.py
```
---

## 🧱 Environment (optional Conda setup)

You can reproduce the environment with:

```yaml
# environment.yml
name: scibd_analysis
channels:
  - conda-forge
dependencies:
  - python=3.9
  - scanpy
  - pandas
  - numpy
  - scipy
```

Activate:

```bash
conda env create -f environment.yml
conda activate scibd_analysis
```
