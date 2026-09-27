#!/usr/bin/env python3
"""Backward-compatible entry point for configurable gene plotting."""

from scrna_expression_eval.plot_cli import main


if __name__ == "__main__":
    raise SystemExit(main())
