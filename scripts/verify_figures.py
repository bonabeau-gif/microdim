#!/usr/bin/env python3
"""Confirm that all data-driven figures for the manuscript are rendered."""
from pathlib import Path
root=Path(__file__).resolve().parents[1]
keys=["Figure2","Figure3","Figure4","Figure5","Supplementary_Figure_S1","Supplementary_Figure_S2","Supplementary_Figure_S3"]
missing=[str(root/'figures'/'rebuilt'/f'{n}_rebuilt.{ext}') for n in keys for ext in ("pdf","png") if not (root/'figures'/'rebuilt'/f'{n}_rebuilt.{ext}').is_file()]
if missing: raise SystemExit('Missing figure(s): '+', '.join(missing))
print('PASS: 7 data-driven figures, each as PNG and PDF')
