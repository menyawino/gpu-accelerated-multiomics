# Processed input layout for reproducible reproduction/validation pipeline

The processed-data pipeline expects tab-delimited files under `data/processed/`.

## Required structure

- `data/processed/gse123976/rna_expression_matrix.tsv`
  - Gene-by-sample matrix (`gene_id` column + sample columns; counts/TPM accepted).
- `data/processed/gse123976/rna_metadata.tsv`
  - Required columns: `sample_id`, `group` (labels include `HF`, `NF` by default).
- `data/processed/gse123976/methylation_gene_table.tsv`
  - Required columns (default names): `gene_id`, `delta_methylation`, `pvalue`.
  - Can be precomputed CpG/DMR-to-gene promoter summary from processed WGBS tables.

- `data/processed/gse197670/beta_matrix.tsv`
  - Probe-by-sample beta matrix (`probe_id` + sample columns).
- `data/processed/gse197670/sample_metadata.tsv`
  - Required columns: `sample_id`, `group` (`Control`, `ICM`, `NICM` by default).
- `data/processed/gse197670/probe_annotation.tsv`
  - Required columns: `probe_id`, `gene_id`, `is_promoter` (TRUE/FALSE).

- `data/processed/gse197672/`
  - Optional companion cohort in same format as GSE197670 for harmonized preprocessing/reference.

- `data/processed/gse116250/rna_expression_matrix.tsv`
  - Gene-by-sample matrix (`gene_id` + sample columns).
- `data/processed/gse116250/sample_metadata.tsv`
  - Required columns: `sample_id`, `group` (`NF`, `DCM`, `ICM`).

- `data/processed/reference/pathways.gmt`
  - GMT pathways used for enrichment.

## Downloading processed matrices from GEO

Use helper script when direct supplementary URLs are available:

```bash
python scripts/download_processed_geo.py --dataset gse123976 --outdir data/processed/gse123976 --url <GEO_SUPPL_URL>
python scripts/download_processed_geo.py --dataset gse197670 --outdir data/processed/gse197670 --url <GEO_SUPPL_URL>
python scripts/download_processed_geo.py --dataset gse116250 --outdir data/processed/gse116250 --url <GEO_SUPPL_URL>
```

If GEO requires manual download, place files in the paths above and rename to the expected filenames.

## Notes

- Validation is at gene/program/pathway level (array vs WGBS platforms differ).
- `GSE197671` is intentionally excluded from myocardial validation (iPSC perturbation RNA-seq).
