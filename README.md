# scRNA-seq Gene Expression Evaluation

This repository provides a reproducible, sample-level workflow for exploring disease-associated gene-expression differences across major cell types in single-cell colon data.

The maintained analysis intentionally separates **cells** from **independent biological samples**. Cells from the same donor/sample are not treated as independent replicates for inferential statistics.

> **Important:** the workflow is exploratory research software. It does not establish causal disease genes, treatment response, or clinical biomarkers.

## The main methodological upgrade

The historical implementation ran Kruskal-Wallis and Mann-Whitney tests directly across individual cells. In single-cell data, many cells can come from the same donor, so treating every cell as an independent replicate can produce **pseudoreplication** and overly small p-values.

The upgraded workflow requires an `adata.obs` column identifying the biological sample/donor. For each cell type, it first computes a **sample-level mean expression vector** across that sample's cells. Statistical tests are then performed across these sample-level aggregates.

```text
cells
  |
  +-- group by biological sample + cell type
  |
  v
sample-level mean expression
  |
  +-- disease groups
  |
  +-- Kruskal-Wallis across samples
  |
  +-- pairwise Mann-Whitney across samples
  |
  +-- BH FDR correction
  |
  v
cell-type-specific results
```

This is closer to the correct replication unit than cell-level testing. It is still not a full count-based pseudobulk differential-expression model such as edgeR/DESeq2/limma-voom.

## Research question

> Which genes show disease-associated differences in sample-level cell-type expression, and which genes separate a target disease group from all other sufficiently represented disease groups?

## Input

Default input:

```text
data/input.h5ad
```

Required `adata.obs` columns are configurable, with defaults:

- `sample_id` — independent biological sample/donor;
- `disease`;
- `major_cluster`.

The analysis rejects a sample ID that maps to more than one disease.

Expression can come from:

- `adata.X` by default;
- a named `adata.layers[...]` via `--layer`;
- `adata.raw.X` via `--use-raw`.

The workflow does not silently transform or normalize the matrix. The biological meaning of the effect-size threshold therefore depends on the expression scale stored in the chosen AnnData matrix.

See [DATA.md](DATA.md).

## Statistical workflow

For each major cell type:

1. aggregate cells to a sample-level mean expression vector;
2. record per-disease descriptive statistics across biological samples;
3. run Kruskal-Wallis tests across disease groups that meet the minimum sample count;
4. apply Benjamini-Hochberg FDR across genes within that cell type;
5. compute a descriptive effect range:
   [
   max(	ext{disease-group mean}) - min(	ext{disease-group mean})
   ]
6. run two-sided Mann-Whitney U tests for each disease-pair contrast;
7. apply Benjamini-Hochberg FDR across genes separately within each disease-pair contrast and cell type;
8. optionally identify genes for which the target disease passes the omnibus FDR/effect criteria and **every eligible target-vs-other pairwise contrast** passes pairwise FDR.

The hard-coded historical denominator `161,060` is no longer used for significance testing. Multiple-testing correction is derived from the tests actually run.

See [METHODS.md](METHODS.md).

## Installation

```bash
git clone https://github.com/seirana/scRNAseq_gene_expression_evaluation.git
cd scRNAseq_gene_expression_evaluation

python -m venv .venv
source .venv/bin/activate

python -m pip install --upgrade pip
python -m pip install -e .
```

For development:

```bash
python -m pip install -e ".[dev]"
```

The maintained workflow uses `anndata` directly rather than requiring the full Scanpy package just to read `.h5ad` files.

## Run the analysis

Example:

```bash
scrna-expression-evaluate \
  --input data/input.h5ad \
  --sample-col sample_id \
  --disease-col disease \
  --cluster-col major_cluster \
  --target-disease UC_inflamed \
  --min-samples-per-group 3 \
  --fdr-threshold 0.05 \
  --effect-threshold 0.25 \
  --outdir outputs
```

The historical entry point remains available:

```bash
python main.py ...
```

If the dataset uses another donor/sample column, pass it explicitly:

```bash
scrna-expression-evaluate \
  --input data/input.h5ad \
  --sample-col donor_id
```

## Outputs

```text
outputs/
├── cell_count_per_disease_and_celltype.csv
├── sample_count_per_disease.csv
├── cell_type_file_map.csv
├── analysis_summary.csv
├── run_metadata.json
├── pseudobulk/
│   ├── <celltype>_sample_means.csv.gz
│   └── <celltype>_analysis_metadata.json
├── stats_per_celltype/
│   └── <celltype>_stats.csv.gz
├── kw_per_celltype/
│   └── <celltype>_KW.csv.gz
├── effectsize_per_celltype/
│   └── <celltype>_effectsize.csv.gz
├── posthoc_per_celltype/
│   └── <celltype>_pairwise_posthoc.csv.gz
└── target_separated/
    └── <celltype>_<target>_separated.csv
```

Each cell-type metadata file records which disease groups met the minimum sample threshold and the FDR correction scope.

`run_metadata.json` records the command, software versions, selected expression source, Git commit when available, and a SHA-256 checksum of the input `.h5ad` file.

## Gene-specific plots

The old plotting script was hard-coded to a local `/home/.../scIBD_Colon.h5ad` path and to `SPOCD1`. It is now configurable:

```bash
scrna-expression-plot-gene \
  --input data/input.h5ad \
  --sample-col sample_id \
  --gene SPOCD1 \
  --outdir plots
```

The plots use **sample-level cell-type means**, not individual cells, and produce:

- one box/jitter plot per cell type;
- a disease × cell-type heatmap;
- a descriptive effect-range bar plot;
- the underlying sample-level values as CSV.

The historical `python plots.py ...` command remains as a wrapper.

## Testing

```bash
python -m pytest
python -m ruff check src tests main.py plots.py
```

Tests cover:

- sample-to-disease consistency;
- Benjamini-Hochberg correction;
- dense and sparse sample-level aggregation;
- sample-level Kruskal-Wallis and Mann-Whitney inputs;
- target-disease selection requiring every eligible contrast;
- end-to-end `.h5ad` analysis;
- configurable gene plotting.

GitHub Actions runs tests and linting on Python 3.10, 3.11, and 3.12 and builds the Docker image.

## Docker

```bash
docker build -t scrna-expression-evaluation .
mkdir -p outputs

docker run --rm \
  -v "$PWD/data:/app/data:ro" \
  -v "$PWD/outputs:/app/outputs" \
  scrna-expression-evaluation \
  --input /app/data/input.h5ad \
  --sample-col sample_id \
  --outdir /app/outputs
```

## Historical results

The repository originally contained many generated CSV and PNG outputs under `results/`. Those run-specific files are not part of the maintained source tree and are excluded from future version control.

Historical PDF material is retained for project context, but the maintained code and current documentation define the upgraded methodology.

## Interpretation limits

A significant result here means that sample-level aggregated expression differs under the specified non-parametric analysis and multiple-testing procedure.

It does **not** by itself establish:

- cell-intrinsic causality;
- a validated disease biomarker;
- treatment response;
- independence from donor covariates, batch, medication, sex, age, or sequencing depth;
- replication in another cohort.

For publication-scale differential expression from raw counts, a donor-aware count model or a suitable mixed model should be considered.

## License

No explicit license file is currently included. Repository visibility alone does not grant reuse rights.
