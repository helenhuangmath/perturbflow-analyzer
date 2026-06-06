# Example Data

`adamson_2016_upr_360x1000.h5ad` is a small real-data subset for smoke-testing
PerturbFlow on a laptop or login node.

Source:

- scPerturb Zenodo record `13350497`
- source file: `AdamsonWeissman2016_GSM2406675_10X001.h5ad`
- source MD5: `232f7e3756d41602bbe434b50662a76f`
- original study: Adamson et al., Cell 2016

The subset keeps 60 cells from each of six labels: one non-targeting/control-like
guide barcode and five UPR perturbations. It keeps the top 1,000 detected genes
and adds a `guide_gene` observation column for the quick-start `prepare` command.

Rebuild it from the public source:

```bash
python scripts/make_example_subset.py
```

The script caches the full downloaded source file under the system temporary
directory unless `--cache-dir` is provided.
