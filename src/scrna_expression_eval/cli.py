from __future__ import annotations

import argparse
import hashlib
import json
import platform
import re
import subprocess
import sys
from importlib.metadata import version
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy

from .core import (
    AnalysisConfig,
    build_pseudobulk_means,
    effect_size_table,
    kruskal_wallis_table,
    pairwise_mannwhitney_table,
    select_target_separated_genes,
    summarize_expression,
    validate_metadata,
)


def _save_json(
    payload: dict[str, object],
    path: Path,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(
        path.suffix + ".tmp"
    )
    temporary.write_text(
        json.dumps(
            payload,
            indent=2,
            sort_keys=True,
        ),
        encoding="utf-8",
    )
    temporary.replace(path)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(
            lambda: handle.read(1024 * 1024),
            b"",
        ):
            digest.update(chunk)
    return digest.hexdigest()


def _git_commit() -> str | None:
    try:
        result = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            check=True,
            capture_output=True,
            text=True,
        )
    except (
        FileNotFoundError,
        subprocess.CalledProcessError,
    ):
        return None
    value = result.stdout.strip()
    return value or None


def _slug(value: str) -> str:
    normalized = re.sub(
        r"[^A-Za-z0-9._-]+",
        "_",
        str(value).strip(),
    ).strip("_")
    return normalized or "unnamed"


def _expression_source(
    adata: ad.AnnData,
    *,
    layer: str | None,
    use_raw: bool,
) -> tuple[object, list[str], str]:
    if layer and use_raw:
        raise ValueError(
            "--layer and --use-raw are mutually exclusive."
        )

    if use_raw:
        if adata.raw is None:
            raise ValueError(
                "--use-raw requested but AnnData.raw is not available."
            )
        return (
            adata.raw.X,
            [str(name) for name in adata.raw.var_names],
            "raw.X",
        )

    if layer:
        if layer not in adata.layers:
            raise ValueError(
                f"Layer {layer!r} not found. "
                f"Available layers: {list(adata.layers.keys())}"
            )
        return (
            adata.layers[layer],
            [str(name) for name in adata.var_names],
            f"layers[{layer!r}]",
        )

    return (
        adata.X,
        [str(name) for name in adata.var_names],
        "X",
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Evaluate disease-associated single-cell expression "
            "at the biological-sample level using per-cell-type "
            "pseudobulk means."
        )
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=Path("data/input.h5ad"),
        help="Input AnnData .h5ad file.",
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=Path("outputs"),
    )
    parser.add_argument(
        "--sample-col",
        default="sample_id",
        help=(
            "obs column identifying independent biological samples/donors. "
            "Inference is performed across these units, not individual cells."
        ),
    )
    parser.add_argument(
        "--disease-col",
        default="disease",
    )
    parser.add_argument(
        "--cluster-col",
        default="major_cluster",
    )
    parser.add_argument(
        "--target-disease",
        default="UC_inflamed",
    )
    parser.add_argument(
        "--min-samples-per-group",
        type=int,
        default=3,
    )
    parser.add_argument(
        "--fdr-threshold",
        type=float,
        default=0.05,
    )
    parser.add_argument(
        "--effect-threshold",
        type=float,
        default=0.25,
        help=(
            "Minimum range of disease-group pseudobulk means. "
            "Its biological meaning depends on the scale stored in AnnData."
        ),
    )
    parser.add_argument(
        "--layer",
        default=None,
        help="Optional AnnData layer to analyze instead of X.",
    )
    parser.add_argument(
        "--use-raw",
        action="store_true",
        help="Use AnnData.raw.X instead of X.",
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()
    config = AnalysisConfig(
        sample_col=args.sample_col,
        disease_col=args.disease_col,
        cluster_col=args.cluster_col,
        target_disease=args.target_disease,
        min_samples_per_group=args.min_samples_per_group,
        fdr_threshold=args.fdr_threshold,
        effect_threshold=args.effect_threshold,
    )

    if not args.input.exists():
        raise FileNotFoundError(
            f"Input AnnData file not found: {args.input}"
        )

    args.outdir.mkdir(
        parents=True,
        exist_ok=True,
    )
    directories = {
        name: args.outdir / name
        for name in [
            "pseudobulk",
            "stats_per_celltype",
            "kw_per_celltype",
            "effectsize_per_celltype",
            "posthoc_per_celltype",
            "target_separated",
        ]
    }
    for directory in directories.values():
        directory.mkdir(
            parents=True,
            exist_ok=True,
        )

    print(
        f"Loading AnnData from: {args.input}"
    )
    adata = ad.read_h5ad(args.input)
    metadata = validate_metadata(
        adata.obs,
        sample_col=config.sample_col,
        disease_col=config.disease_col,
        cluster_col=config.cluster_col,
    )
    matrix, gene_names, source_name = (
        _expression_source(
            adata,
            layer=args.layer,
            use_raw=args.use_raw,
        )
    )

    if len(set(gene_names)) != len(gene_names):
        raise ValueError(
            "Gene names must be unique for reproducible tabular outputs."
        )

    cell_counts = (
        metadata[
            [
                config.cluster_col,
                config.disease_col,
            ]
        ]
        .value_counts()
        .reset_index(name="n_cells")
        .sort_values(
            [
                config.cluster_col,
                config.disease_col,
            ],
            kind="stable",
        )
    )
    cell_counts.to_csv(
        args.outdir
        / "cell_count_per_disease_and_celltype.csv",
        index=False,
    )

    sample_disease = (
        metadata[
            [
                config.sample_col,
                config.disease_col,
            ]
        ]
        .drop_duplicates()
    )
    (
        sample_disease[
            config.disease_col
        ]
        .value_counts()
        .sort_index()
        .rename("n_biological_samples")
        .reset_index()
        .to_csv(
            args.outdir
            / "sample_count_per_disease.csv",
            index=False,
        )
    )

    cell_types = sorted(
        metadata[
            config.cluster_col
        ].unique().tolist()
    )

    file_map_rows: list[
        dict[str, str]
    ] = []
    summary_rows: list[
        dict[str, object]
    ] = []

    for cell_type in cell_types:
        print(
            f"Processing cell type: {cell_type}"
        )
        slug = _slug(cell_type)
        file_map_rows.append(
            {
                "cell_type": cell_type,
                "file_slug": slug,
            }
        )

        pseudobulk = build_pseudobulk_means(
            matrix,
            metadata,
            gene_names,
            sample_col=config.sample_col,
            disease_col=config.disease_col,
            cluster_col=config.cluster_col,
            cluster_value=cell_type,
        )

        pseudobulk_frame = pd.DataFrame(
            pseudobulk.expression,
            columns=pseudobulk.gene_names,
        )
        pseudobulk_frame.insert(
            0,
            "n_cells",
            pseudobulk.metadata[
                "n_cells"
            ].to_numpy(),
        )
        pseudobulk_frame.insert(
            0,
            config.disease_col,
            pseudobulk.metadata[
                config.disease_col
            ].to_numpy(),
        )
        pseudobulk_frame.insert(
            0,
            config.sample_col,
            pseudobulk.metadata[
                config.sample_col
            ].to_numpy(),
        )
        pseudobulk_frame.to_csv(
            directories["pseudobulk"]
            / f"{slug}_sample_means.csv.gz",
            index=False,
            compression="gzip",
        )

        stats = summarize_expression(
            pseudobulk,
            disease_col=config.disease_col,
        )
        stats.to_csv(
            directories[
                "stats_per_celltype"
            ]
            / f"{slug}_stats.csv.gz",
            index=False,
            compression="gzip",
        )

        try:
            kw, kw_metadata = (
                kruskal_wallis_table(
                    pseudobulk,
                    disease_col=(
                        config.disease_col
                    ),
                    min_samples_per_group=(
                        config.min_samples_per_group
                    ),
                )
            )
            effects = effect_size_table(
                pseudobulk,
                disease_col=config.disease_col,
                min_samples_per_group=(
                    config.min_samples_per_group
                ),
            )
            pairwise, pairwise_metadata = (
                pairwise_mannwhitney_table(
                    pseudobulk,
                    disease_col=(
                        config.disease_col
                    ),
                    min_samples_per_group=(
                        config.min_samples_per_group
                    ),
                )
            )
        except ValueError as error:
            summary_rows.append(
                {
                    "cell_type": cell_type,
                    "n_cells": int(
                        pseudobulk.metadata[
                            "n_cells"
                        ].sum()
                    ),
                    "n_biological_samples": int(
                        len(
                            pseudobulk.metadata
                        )
                    ),
                    "n_eligible_diseases": 0,
                    "n_kw_fdr_significant": 0,
                    "n_target_separated": 0,
                    "status": str(error),
                }
            )
            continue

        kw.to_csv(
            directories["kw_per_celltype"]
            / f"{slug}_KW.csv.gz",
            index=False,
            compression="gzip",
        )
        effects.to_csv(
            directories[
                "effectsize_per_celltype"
            ]
            / f"{slug}_effectsize.csv.gz",
            index=False,
            compression="gzip",
        )
        pairwise.to_csv(
            directories[
                "posthoc_per_celltype"
            ]
            / f"{slug}_pairwise_posthoc.csv.gz",
            index=False,
            compression="gzip",
        )

        selected = (
            select_target_separated_genes(
                kw,
                effects,
                pairwise,
                target_disease=(
                    config.target_disease
                ),
                eligible_diseases=(
                    kw_metadata[
                        "eligible_diseases"
                    ]
                ),
                fdr_threshold=(
                    config.fdr_threshold
                ),
                effect_threshold=(
                    config.effect_threshold
                ),
            )
        )
        selected.to_csv(
            directories[
                "target_separated"
            ]
            / (
                f"{slug}_"
                f"{_slug(config.target_disease)}"
                "_separated.csv"
            ),
            index=False,
        )

        _save_json(
            {
                "cell_type": cell_type,
                "kw": kw_metadata,
                "pairwise": (
                    pairwise_metadata
                ),
            },
            directories["pseudobulk"]
            / f"{slug}_analysis_metadata.json",
        )

        summary_rows.append(
            {
                "cell_type": cell_type,
                "n_cells": int(
                    pseudobulk.metadata[
                        "n_cells"
                    ].sum()
                ),
                "n_biological_samples": int(
                    len(
                        pseudobulk.metadata
                    )
                ),
                "n_eligible_diseases": int(
                    len(
                        kw_metadata[
                            "eligible_diseases"
                        ]
                    )
                ),
                "n_kw_fdr_significant": int(
                    (
                        kw["FDR"]
                        < config.fdr_threshold
                    ).sum()
                ),
                "n_target_separated": int(
                    len(selected)
                ),
                "status": "ok",
            }
        )

    pd.DataFrame(
        file_map_rows
    ).to_csv(
        args.outdir
        / "cell_type_file_map.csv",
        index=False,
    )
    summary = pd.DataFrame(
        summary_rows
    )
    summary.to_csv(
        args.outdir
        / "analysis_summary.csv",
        index=False,
    )

    run_metadata = {
        "command": (
            "scrna-expression-evaluate"
        ),
        "arguments": {
            key: (
                str(value)
                if isinstance(value, Path)
                else value
            )
            for key, value
            in vars(args).items()
        },
        "method": {
            "inference_unit": (
                "biological sample/donor"
            ),
            "aggregation": (
                "mean expression across cells "
                "within sample and cell type"
            ),
            "omnibus_test": (
                "Kruskal-Wallis across "
                "eligible disease groups"
            ),
            "omnibus_fdr_scope": (
                "genes within cell type"
            ),
            "pairwise_test": (
                "two-sided Mann-Whitney U"
            ),
            "pairwise_fdr_scope": (
                "genes within each disease "
                "contrast and cell type"
            ),
            "target_selection": {
                "target_disease": (
                    config.target_disease
                ),
                "kw_fdr_lt": (
                    config.fdr_threshold
                ),
                "effect_size_gt": (
                    config.effect_threshold
                ),
                "all_target_pairwise_fdr_lt": (
                    config.fdr_threshold
                ),
                "comparison_population": (
                    "disease groups meeting "
                    "min_samples_per_group"
                ),
            },
        },
        "input": {
            "path": str(
                args.input.resolve()
            ),
            "sha256": _sha256(
                args.input
            ),
            "expression_source": (
                source_name
            ),
            "n_cells": int(
                adata.n_obs
            ),
            "n_genes": int(
                len(gene_names)
            ),
            "n_biological_samples": int(
                metadata[
                    config.sample_col
                ].nunique()
            ),
        },
        "software": {
            "python": (
                sys.version.split()[0]
            ),
            "platform": (
                platform.platform()
            ),
            "anndata": version("anndata"),
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "scipy": scipy.__version__,
        },
        "git_commit": _git_commit(),
    }
    _save_json(
        run_metadata,
        args.outdir
        / "run_metadata.json",
    )

    print(
        "Done. Results written to: "
        f"{args.outdir.resolve()}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
