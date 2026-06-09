# Outputs

Each run writes a result folder with reports, tables, plots, a final AnnData file, and a bundle for web viewers.

```text
results/my_run/
├── interactive_report.html
├── report.html
├── final_adata.h5ad
├── checkpoint.json
├── summary.json
├── csv/
├── plots/
├── bundle/
└── agent_handoff/
```

## Important Files

| Path | Purpose |
| --- | --- |
| `interactive_report.html` | Browser-based report for scientific review |
| `csv/deg_summary.csv` | One-row-per-perturbation DEG summary |
| `csv/cscore_summary.csv` | Connectivity rewiring scores |
| `csv/effect_decomposition.csv` | Transcriptional and state-shift effect summaries |
| `bundle/manifest.json` | Versioned artifact index for downstream viewers |
| `agent_handoff/interpretation_context.md` | Agent interpretation context |
