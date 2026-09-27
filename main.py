#!/usr/bin/env python3
"""Backward-compatible entry point for sample-level scRNA expression evaluation."""

from scrna_expression_eval.cli import main


if __name__ == "__main__":
    raise SystemExit(main())
