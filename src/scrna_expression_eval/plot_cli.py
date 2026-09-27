from __future__ import annotations

import argparse
import re
from pathlib import Path

import anndata as ad
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .core import (
    build_pseudobulk_means,
    validate_metadata,
)


def _slug(value: str) -> str:
    result = re.sub(
        r"[^A-Za-z0-9._-]+",
        "_",
        str(value).strip(),
    ).strip("_")
    return result or "unnamed"


def _expression_source(
    adata: ad.AnnData,
    layer: str | None,
    use_raw: bool,
) -> tuple[object, list[str]]:
    if layer and use_raw:
        raise ValueError(
            "--layer and --use-raw are mutually exclusive."
        )
    if use_raw:
        if adata.raw is None:
            raise ValueError(
                "--use-raw requested but AnnData.raw is unavailable."
            )
        return (
            adata.raw.X,
            [str(value) for value in adata.raw.var_names],
        )
    if layer:
        if layer not in adata.layers:
            raise ValueError(
                f"Layer {layer!r} not found."
            )
        return (
            adata.layers[layer],
            [str(value) for value in adata.var_names],
        )
    return (
        adata.X,
        [str(value) for value in adata.var_names],
    )


def _ordered_diseases(
    observed: list[str],
    requested: list[str],
) -> list[str]:
    ordered = [
        disease
        for disease in requested
        if disease in observed
    ]
    ordered.extend(
        disease
        for disease in sorted(observed)
        if disease not in ordered
    )
    return ordered


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Plot one gene using sample-level "
            "per-cell-type pseudobulk means."
        )
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=Path("data/input.h5ad"),
    )
    parser.add_argument(
        "--gene",
        default="SPOCD1",
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=Path("plots"),
    )
    parser.add_argument(
        "--sample-col",
        default="sample_id",
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
        "--disease-order",
        default=(
            "Healthy,UC_non_inflamed,"
            "UC_inflamed,CD_inflamed,"
            "Colitis_inflamed"
        ),
    )
    parser.add_argument(
        "--layer",
        default=None,
    )
    parser.add_argument(
        "--use-raw",
        action="store_true",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()
    if not args.input.exists():
        raise FileNotFoundError(
            f"Input AnnData file not found: {args.input}"
        )

    args.outdir.mkdir(
        parents=True,
        exist_ok=True,
    )
    adata = ad.read_h5ad(
        args.input
    )
    metadata = validate_metadata(
        adata.obs,
        sample_col=args.sample_col,
        disease_col=args.disease_col,
        cluster_col=args.cluster_col,
    )
    matrix, gene_names = _expression_source(
        adata,
        args.layer,
        args.use_raw,
    )

    if args.gene not in gene_names:
        raise ValueError(
            f"Gene {args.gene!r} not found in selected expression source."
        )
    gene_index = gene_names.index(
        args.gene
    )

    requested_order = [
        item.strip()
        for item in args.disease_order.split(",")
        if item.strip()
    ]
    rng = np.random.default_rng(
        args.seed
    )

    records: list[
        dict[str, object]
    ] = []

    for cell_type in sorted(
        metadata[
            args.cluster_col
        ].unique().tolist()
    ):
        pseudobulk = build_pseudobulk_means(
            matrix,
            metadata,
            gene_names,
            sample_col=args.sample_col,
            disease_col=args.disease_col,
            cluster_col=args.cluster_col,
            cluster_value=cell_type,
        )
        values = pseudobulk.expression[
            :,
            gene_index,
        ]
        frame = pseudobulk.metadata[
            [
                args.sample_col,
                args.disease_col,
                "n_cells",
            ]
        ].copy()
        frame["mean_expression"] = (
            values
        )
        frame["cell_type"] = (
            cell_type
        )
        records.extend(
            frame.to_dict(
                orient="records"
            )
        )

        observed = (
            frame[
                args.disease_col
            ]
            .astype(str)
            .unique()
            .tolist()
        )
        order = _ordered_diseases(
            observed,
            requested_order,
        )
        grouped = [
            frame.loc[
                frame[
                    args.disease_col
                ].astype(str)
                == disease,
                "mean_expression",
            ].to_numpy(dtype=float)
            for disease in order
        ]

        plt.figure(
            figsize=(
                max(7, len(order) * 1.2),
                4.5,
            )
        )
        plt.boxplot(
            grouped,
            tick_labels=order,
            showfliers=False,
        )
        for position, values_group in enumerate(
            grouped,
            start=1,
        ):
            jitter = rng.normal(
                position,
                0.04,
                size=len(values_group),
            )
            plt.scatter(
                jitter,
                values_group,
                alpha=0.7,
                s=18,
            )
        plt.title(
            f"{args.gene} sample-level mean expression\n"
            f"Cell type: {cell_type}"
        )
        plt.xlabel("Disease")
        plt.ylabel(
            "Mean expression within "
            "sample/cell type"
        )
        plt.xticks(
            rotation=30,
            ha="right",
        )
        plt.tight_layout()
        plt.savefig(
            args.outdir
            / (
                f"{_slug(args.gene)}_"
                f"{_slug(cell_type)}"
                "_pseudobulk_boxplot.png"
            ),
            dpi=150,
        )
        plt.close()

    values_frame = pd.DataFrame(
        records
    )
    values_frame.to_csv(
        args.outdir
        / (
            f"{_slug(args.gene)}"
            "_pseudobulk_values.csv"
        ),
        index=False,
    )

    grouped_means = (
        values_frame.groupby(
            [
                "cell_type",
                args.disease_col,
            ],
            sort=True,
        )["mean_expression"]
        .mean()
        .reset_index()
    )
    heat = grouped_means.pivot(
        index="cell_type",
        columns=args.disease_col,
        values="mean_expression",
    )
    ordered_columns = _ordered_diseases(
        heat.columns.astype(str).tolist(),
        requested_order,
    )
    heat = heat.reindex(
        columns=ordered_columns
    )

    plt.figure(
        figsize=(
            max(7, len(heat.columns) * 1.2),
            max(5, len(heat.index) * 0.5),
        )
    )
    image = plt.imshow(
        heat.to_numpy(dtype=float),
        aspect="auto",
    )
    plt.colorbar(
        image,
        label=(
            "Mean of sample-level "
            "mean expression"
        ),
    )
    plt.xticks(
        range(len(heat.columns)),
        heat.columns,
        rotation=30,
        ha="right",
    )
    plt.yticks(
        range(len(heat.index)),
        heat.index,
    )
    plt.title(
        f"{args.gene} across disease and cell type"
    )
    plt.tight_layout()
    plt.savefig(
        args.outdir
        / (
            f"{_slug(args.gene)}"
            "_pseudobulk_heatmap.png"
        ),
        dpi=150,
    )
    plt.close()

    effect = (
        grouped_means.groupby(
            "cell_type"
        )["mean_expression"]
        .agg(
            lambda values: (
                values.max()
                - values.min()
            )
        )
        .sort_values(
            ascending=False
        )
    )
    plt.figure(
        figsize=(
            max(
                8,
                len(effect) * 0.8,
            ),
            4.5,
        )
    )
    plt.bar(
        effect.index,
        effect.values,
    )
    plt.ylabel(
        "Range of disease-group "
        "sample means"
    )
    plt.title(
        f"{args.gene} descriptive effect range"
    )
    plt.xticks(
        rotation=35,
        ha="right",
    )
    plt.tight_layout()
    plt.savefig(
        args.outdir
        / (
            f"{_slug(args.gene)}"
            "_effect_range.png"
        ),
        dpi=150,
    )
    plt.close()

    print(
        "Plots written to: "
        f"{args.outdir.resolve()}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
