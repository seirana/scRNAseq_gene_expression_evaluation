# Historical and generated results

This directory previously contained many generated CSV and PNG files from the historical cell-level workflow.

Those run-specific outputs are removed from the maintained source tree because:

- they can be regenerated from the appropriate source data;
- they substantially increase repository size;
- they were produced with the historical cell-level inferential method;
- keeping generated tables beside source code makes the current methodology ambiguous.

The historical `Methodological Overview.pdf` is retained for project context.

For current methodology, use:

- `../README.md`;
- `../METHODS.md`;
- `../DATA.md`.

New runs should write to `outputs/` or another user-selected output directory and should not be committed unless there is a deliberate archival reason.
