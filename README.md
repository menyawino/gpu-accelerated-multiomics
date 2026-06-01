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

Incremental processing of only currently complete FASTQ pairs, without downstream integrative outputs:

```bash
./run_analysis.sh --available-only --cores 8 --jobs 4
```

Keep watching the raw FASTQ directory and process newly completed pairs as they arrive:

```bash
./run_analysis.sh --available-only --watch-interval 120 --cores 8 --jobs 4
```

This mode always skips the matrix and network build rules. It trims any complete pairs immediately, adds Parabricks surrogate BAM targets when Docker and references are available, and only adds RNA-seq or WGBS per-sample outputs when the required references and run metadata are present.

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

## Reproduce and validate GSE123976 methylation–expression program

This processed-data pipeline reproduces discovery signals from GSE123976 and validates them externally using methylation-array and transcriptome cohorts.

### 1) Prepare processed inputs

Place processed matrices/metadata in the standardized layout described in:

- `data/processed/README.md`

Optional helper for direct GEO supplementary downloads:

```bash
python scripts/download_processed_geo.py --dataset gse123976 --outdir data/processed/gse123976 --url <GEO_SUPPL_URL>
```

### 2) Configure

Pipeline config is in:

- `config/processed_analysis.yaml`

You can edit labels, thresholds, and file paths there.

### 3) Create environment

```bash
conda env create -f environment.yml
conda activate hf_processed
```

### 4) Run end-to-end

```bash
chmod +x run_processed_analysis.sh
./run_processed_analysis.sh --cores 8 --jobs 4
```

Dry run:

```bash
./run_processed_analysis.sh --dry-run
```

### 5) Outputs

- `results/gse123976/`: DEG table, integrated gene classes, enrichment, QC plot
- `results/gse197670/`: DMP table, promoter gene summaries, concordance stats/plot
- `results/gse116250/`: DCM/NF + ICM/NF DEG, GSEA-like enrichment, signature scores/plot
- `results/integration/integrated_gene_summary.tsv`
- `reports/gse123976_validation_report.qmd`

### Workflow + environment

- Workflow: `workflow/processed/Snakefile` (Snakemake)
- Environment: `environment.yml` (also `envs/processed.yaml`)
