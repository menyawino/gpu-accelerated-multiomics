# GPU-Accelerated Multi-Omics Pipeline (WGBS + RNA-seq + Isoform Network)

This project provides a Snakemake workflow to:

1. Discover and download GSE123976 runs from SRA.
2. Process WGBS data to methylation calls.
3. Process RNA-seq data to isoform- and gene-level expression.
4. Build multi-omic correlation networks connecting methylation and expression at:
   - gene level
   - isoform/transcript level

## GPU and Parabricks usage

NVIDIA Parabricks is used in all compatible points through `scripts/pbrun.sh` and the `parabricks_fq2bam_surrogate` rule.

Important technical note:
- There is no native Parabricks bisulfite-aware methylation calling equivalent to Bismark in this workflow.
- For biologically correct WGBS methylation extraction, Bismark is used for alignment + methylation calling.
- RNA isoform quantification is done by STAR + Salmon (standard and biologically appropriate path).

This satisfies a GPU-first architecture while preserving validity for methylation and isoform inference.

## Files

- `setup_tools.sh`: install/check dependencies and GPU stack.
- `run_analysis.sh`: calls setup, then launches Snakemake.
- `workflow/Snakefile`: full pipeline DAG.
- `config/config.yaml`: pipeline configuration.
- `config/manual_pairs.tsv`: optional deterministic WGBS/RNA sample pairing override.

## Prerequisites for T4 execution host

Use a Linux machine with:

1. NVIDIA driver and `nvidia-smi` working.
2. Docker + NVIDIA Container Toolkit (`docker run --gpus all ...`).
3. NGC access (`docker login nvcr.io`) for pulling Parabricks image.

## Configure references

Place references at:

- `resources/ref/GRCh38.primary_assembly.genome.fa`
- `resources/ref/gencode.v44.annotation.gtf`
- `resources/ref/gencode.v44.transcripts.fa`

Adjust paths in `config/config.yaml` if needed.

## Run

```bash
chmod +x setup_tools.sh run_analysis.sh scripts/pbrun.sh scripts/*.py
./run_analysis.sh --cores 32 --jobs 8
```

Dry run:

```bash
./run_analysis.sh --dry-run --cores 8 --jobs 4
```

## Pairing strategy

Automated pairing is inferred from metadata subject-like tokens. For publication-grade analysis, provide explicit mapping in:

- `config/manual_pairs.tsv`

Columns required:

- `subject_id`
- `rnaseq_run`
- `wgbs_run`

## Key outputs

- `results/matrices/gene_expression_tpm.tsv`
- `results/matrices/isoform_expression_tpm.tsv`
- `results/matrices/gene_methylation_beta.tsv`
- `results/matrices/isoform_methylation_beta.tsv`
- `results/network/gene_network_edges.tsv`
- `results/network/isoform_network_edges.tsv`
- `results/network/gene_network.graphml`
- `results/network/isoform_network.graphml`
