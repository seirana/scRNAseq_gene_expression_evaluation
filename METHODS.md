# Methods

## Unit of inference

The maintained workflow treats the **biological sample/donor** as the independent replication unit.

For each cell type and biological sample, expression is aggregated by taking the mean across all cells from that sample in that cell type.

This addresses a major limitation of the historical cell-level testing workflow: cells from the same donor are correlated and should not automatically be counted as independent biological replicates.

## Expression source

The workflow analyzes exactly one AnnData expression source:

- `adata.X`;
- a named `adata.layers[...]`;
- or `adata.raw.X`.

No normalization or transformation is silently applied. Users must know whether the selected values are raw counts, normalized counts, log-transformed expression, or another representation.

Because aggregation uses a mean rather than a raw-count sum, this workflow is best described as a **sample-level cell-type expression analysis**, not a count-model pseudobulk differential-expression pipeline.

## Metadata validation

The following are required:

- a biological-sample column;
- a disease column;
- a cell-type/cluster column.

Missing or empty values are rejected.

Each biological sample must map to exactly one disease. A sample that appears with multiple disease labels causes the run to fail.

## Sample-level aggregation

For one cell type, let (x_{c,g}) be the selected AnnData expression value for cell (c) and gene (g).

For sample (s), the aggregate is:

[
ar{x}_{s,g} = rac{1}{n_s}sum_{c in s} x_{c,g}
]

where (n_s) is the number of cells from sample (s) in that cell type.

The number of contributing cells is retained for QC.

## Minimum sample rule

A disease group is eligible for inferential testing within a cell type only if it contains at least `--min-samples-per-group` biological samples.

The default is 3.

Disease groups below the threshold are retained in descriptive sample-level outputs but excluded from Kruskal-Wallis and pairwise inferential tests. The per-cell-type metadata explicitly records eligible and excluded groups.

## Descriptive statistics

For every gene and disease group, the pipeline reports across biological-sample aggregates:

- number of samples;
- mean;
- sample variance;
- minimum;
- maximum.

## Omnibus test

For each gene within each cell type, a Kruskal-Wallis test compares the eligible disease groups.

If all tested values for a gene are identical, the pipeline records:

- (H = 0);
- (p = 1).

The workflow also reports an epsilon-squared effect statistic:

[
epsilon^2 = maxleft(0, rac{H-k+1}{n-k}ight)
]

where (H) is the Kruskal-Wallis statistic, (k) is the number of disease groups, and (n) is the number of biological samples included in that cell type's test.

## Omnibus multiple-testing correction

Benjamini-Hochberg FDR is applied across all gene-level Kruskal-Wallis p-values within a cell type.

The historical fixed denominator of 161,060 tests is not used.

## Descriptive mean-range effect size

The repository preserves a simple effect measure related to the historical implementation:

[
	ext{effect range} =
max_d(ar{x}_{d,g})
-
min_d(ar{x}_{d,g})
]

where (ar{x}_{d,g}) is the mean of biological-sample aggregates for disease (d).

This value is scale dependent. A threshold such as 0.25 has no universal biological meaning and must be interpreted relative to the expression representation stored in the AnnData object.

## Pairwise post-hoc tests

For each eligible disease pair and gene, the pipeline runs a two-sided Mann-Whitney U test across biological-sample aggregates.

It also reports rank-biserial correlation:

[
r_{rb} = rac{2U}{n_1n_2} - 1
]

using the orientation of disease 1 relative to disease 2.

## Pairwise multiple-testing correction

For each cell type and disease-pair contrast, Benjamini-Hochberg correction is applied across genes.

Therefore each disease pair is treated as its own gene-testing family.

The workflow does not claim formal family-wise error control across every possible disease contrast simultaneously.

## Target-disease separation rule

A gene is placed in the target-separated output when:

1. its Kruskal-Wallis FDR is below `--fdr-threshold`;
2. its descriptive mean-range effect size is above `--effect-threshold`;
3. the target disease is among the eligible disease groups;
4. every eligible target-vs-other disease contrast has pairwise FDR below the same threshold.

This rule applies only to disease groups that meet the minimum biological-sample requirement. It makes no claim about underrepresented disease groups.

## Why this is not a complete differential-expression model

This repository does not currently model:

- raw-count mean/variance relationships;
- library-size offsets;
- donor-level covariates;
- paired/repeated measurements;
- batch effects;
- random effects;
- zero inflation;
- compositional cell-type abundance;
- donor-by-condition interactions.

If raw counts and sufficient metadata are available, publication-scale work should consider established donor-aware pseudobulk count models or statistically appropriate mixed models.

## Multiple cells per sample

Cell counts affect how precisely a sample-level mean is estimated, but the inferential sample size is the number of independent biological samples, not the number of cells.

A donor with 5,000 cells is still one biological replicate.

## Reproducibility

A run records:

- input file path and SHA-256 checksum;
- selected AnnData expression source;
- metadata column names;
- target disease;
- sample-count threshold;
- FDR threshold;
- effect threshold;
- Python/platform information;
- library versions;
- Git commit when available.

Generated result tables are not committed to the repository by default.
