# 🧬 scIBD Colon – Cluster-wise Gene Expression Summary & Disease-specific Significance (Kruskal–Wallis)

This repository provides a reproducible pipeline for analyzing **single-cell colon data (scIBD)** to identify genes whose expression levels differ across diseases and within cell clusters.  

The workflow:
1. Computes per-gene summary statistics for each **disease × minor_cluster**.
2. Aggregates disease-level averages across clusters.
3. Combines all diseases into a unified table.
4. Applies **Kruskal–Wallis tests** to detect significant disease effects on gene expression.
5. Repeats these tests **within each minor_cluster** to account for cell-type-specific differences.

---

## ❓ Research Question

> *Which genes exhibit disease-associated expression differences in single-cell colon data, and how consistent are these effects across minor cell clusters?*

---

## 📥 Input

- **AnnData file:**  
  `/work_beegfs/sukmb655/data_scDRS/sinlge cell datasets/scIBD_Colon.h5ad`

- **Required columns in** `adata.obs`  
  - `disease`
  - `minor_cluster`

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
python run_scibd_pipeline.py
```

This script:

1. Loads the `.h5ad` file using **Scanpy**.
2. Computes per-cluster gene statistics.
3. Aggregates disease-level averages.
4. Merges all diseases into one file.
5. Runs global and per-cluster Kruskal–Wallis tests.
6. Saves all results in organized output folders.

---

## 🧩 Step-by-Step Overview

### **1️⃣ Per-disease × minor_cluster Gene Statistics**

For each disease:

* Subset to cells of that disease.
* For each `minor_cluster`, compute per-gene:

  * **mean**, **variance**, **min**, **max**

➡️ Output:
`cluster_gene_stats/{disease}_gene_stats.csv`

* Rows: genes
* Columns: `{cluster}_mean`, `{cluster}_var`, `{cluster}_min`, `{cluster}_max`

---

### **2️⃣ Per-disease Averages across Clusters**

Compute the **mean of cluster means** for each gene.

➡️ Output:
`disease_avg_expression/{disease}_avg_expression.csv`
Columns: `gene`, `avg_expression_across_clusters`

---

### **3️⃣ Combined Disease Matrix**

Merge all per-disease averages into one file.

➡️ Output:
`disease_avg_expression_all.csv`

* Rows: genes
* Columns: one per disease

---

### **4️⃣ Global Kruskal–Wallis (All Diseases)**

Tests whether each gene’s expression differs between diseases across all cells (ignoring clusters).

➡️ Output:
`gene_disease_significance.csv`
Columns:
| gene | kruskal_p | FDR |

---

### **5️⃣ Per-cluster Kruskal–Wallis**

* Rounds `adata.X` to **3 decimal digits** for numerical stability.
* Runs Kruskal–Wallis tests for each gene **within each cluster**, comparing diseases only inside that cluster.
* Applies **Benjamini–Hochberg FDR correction**.

➡️ Output (in `kruskal_per_cluster/`):

* `{cluster}_gene_disease_significance.csv`

  * Columns: `gene`, `H_statistic`, `p_value`, `FDR`
* `{cluster}_group_sizes.csv`

  * Columns: `disease`, `n_cells`

---

## 🧪 Statistical Methods

| Test                                     | Purpose                                                           | Scope                 |
| ---------------------------------------- | ----------------------------------------------------------------- | --------------------- |
| **Kruskal–Wallis (scipy.stats.kruskal)** | Detects differences in median gene expression between ≥2 diseases | Global or per-cluster |
| **Benjamini–Hochberg (FDR)**             | Controls for multiple testing                                     | Applied per test set  |

> The per-cluster test controls for cluster composition differences between diseases.

---

## 🎯 Precision and Rounding

| Data Type           | Range | Rounding | Example          |
| ------------------- | ----- | -------- | ---------------- |
| Expression values   | 0–10  | 3 digits | 2.456 → 2.456    |
| Kruskal H statistic | 0–100 | 3 digits | 12.349 → 12.349  |
| p-values / FDR      | 0–1   | 3 digits | 0.000234 → 0.000 |

Expression rounding is done safely for both sparse and dense matrices:

```python
if sp.issparse(adata.X):
    X = adata.X.tocoo(copy=True)
    X.data = np.round(X.data, 3)
    adata.X = X.tocsr()
else:
    adata.X = np.round(adata.X, 3)
```

---

## 📂 Output Structure

```
.
├── run_scibd_pipeline.py
├── cluster_gene_stats/
│   ├── UC_inflamed_gene_stats.csv
│   ├── Healthy_gene_stats.csv
│   └── ...
├── disease_avg_expression/
│   ├── UC_inflamed_avg_expression.csv
│   ├── Healthy_avg_expression.csv
│   └── ...
├── disease_avg_expression_all.csv
├── gene_disease_significance.csv
└── kruskal_per_cluster/
    ├── Enterocyte_gene_disease_significance.csv
    ├── Enterocyte_group_sizes.csv
    └── ...
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

```
