# Data requirements and provenance

## Input format

The maintained pipeline accepts an AnnData `.h5ad` file.

Default path:

```text
data/input.h5ad
```

The input data are not committed to this public repository.

## Required observation metadata

The default expected `adata.obs` columns are:

```text
sample_id
disease
major_cluster
```

All three names are configurable from the CLI.

### Biological sample ID

The sample column must identify the independent biological replicate: typically a donor, patient, biopsy, or other genuinely independent sampling unit appropriate to the study design.

A technical barcode, cell ID, or cluster label is **not** an appropriate replacement.

### Disease label

Each biological sample must map to exactly one disease label in the maintained workflow.

If the same donor contributes repeated conditions or timepoints, a more explicit longitudinal/paired statistical design is needed rather than silently treating those observations as independent.

### Cell-type label

The cluster column should contain the cell-type or major-cluster assignment to analyze separately.

## Expression matrix provenance

Record which AnnData matrix is used:

- `X`;
- a named layer;
- or `raw.X`.

Also record how that matrix was created, including where applicable:

- raw count source;
- normalization method;
- size-factor method;
- log transformation;
- gene filtering;
- batch correction;
- integration procedure;
- imputation or smoothing.

The maintained code does not infer these preprocessing choices from the matrix.

## Study-level provenance

For a reproducible analysis, preserve:

- dataset accession/source;
- cohort inclusion/exclusion rules;
- donor/sample identifiers;
- disease definitions;
- tissue and anatomical site;
- sequencing chemistry/platform;
- genome/transcriptome reference;
- gene annotation version;
- preprocessing pipeline/version;
- QC thresholds;
- cell-type annotation method;
- batch-correction/integration method;
- treatment status and major clinical covariates when relevant.

## Privacy

Single-cell datasets can contain sensitive genomic and phenotype information.

Do not commit controlled-access or identifiable donor-level data to the public repository. Keep them in an approved local or controlled environment and mount/copy only the required input at runtime.

## Generated outputs

Outputs can include donor/sample identifiers. Treat generated CSV files according to the same data-governance requirements as the source dataset.

The repository ignores `outputs/`, `plots/`, and future generated `results/` files by default.
