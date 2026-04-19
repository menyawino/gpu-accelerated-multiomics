#!/usr/bin/env python3
import argparse
import shutil
import subprocess
import sys
import time
from pathlib import Path

import yaml


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", required=True, type=Path)
    parser.add_argument("--configfile", required=True, type=Path)
    parser.add_argument("--cores", required=True, type=int)
    parser.add_argument("--jobs", required=True, type=int)
    parser.add_argument("--profile", required=True)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--watch-interval", type=int, default=0)
    return parser.parse_args()


def load_config(configfile: Path) -> dict:
    with configfile.open() as handle:
        return yaml.safe_load(handle)


def resolve_path(project_root: Path, raw_path: str) -> Path:
    path = Path(raw_path)
    return path if path.is_absolute() else project_root / path


def nonempty_file(path: Path) -> bool:
    return path.is_file() and path.stat().st_size > 0


def discover_complete_runs(raw_dir: Path) -> tuple[list[str], list[str]]:
    complete = []
    incomplete = []

    for r1 in sorted(raw_dir.glob("*_1.fastq.gz")):
        run = r1.name[: -len("_1.fastq.gz")]
        r2 = raw_dir / f"{run}_2.fastq.gz"
        if nonempty_file(r1) and nonempty_file(r2):
            complete.append(run)
        else:
            incomplete.append(run)

    return complete, sorted(set(incomplete))


def maybe_discover_metadata(project_root: Path, config: dict, metadata_dir: Path) -> None:
    rnaseq_list = metadata_dir / "rnaseq_runs.txt"
    wgbs_list = metadata_dir / "wgbs_runs.txt"
    if rnaseq_list.exists() and wgbs_list.exists():
        return

    metadata_dir.mkdir(parents=True, exist_ok=True)
    manual_pairs = resolve_path(project_root, config["pairing"]["manual_pairs_tsv"])
    cmd = [
        sys.executable,
        str(project_root / "scripts" / "discover_gse123976_runs.py"),
        "--gse",
        str(config["project"]["gse_accession"]),
        "--metadata-dir",
        str(metadata_dir),
        "--manual-pairs",
        str(manual_pairs),
    ]
    try:
        subprocess.run(cmd, cwd=project_root, check=True)
    except subprocess.CalledProcessError as exc:
        print(f"Metadata discovery failed: {exc}", file=sys.stderr)


def read_run_list(path: Path) -> set[str]:
    if not path.exists():
        return set()
    return {line.strip() for line in path.read_text().splitlines() if line.strip()}


def build_targets(project_root: Path, config: dict, complete_runs: list[str]) -> tuple[list[str], list[str]]:
    warnings = []
    outdir = config["resources"]["outdir"]
    metadata_dir = resolve_path(project_root, config["resources"]["metadata_dir"])
    maybe_discover_metadata(project_root, config, metadata_dir)

    rnaseq_runs = read_run_list(metadata_dir / "rnaseq_runs.txt")
    wgbs_runs = read_run_list(metadata_dir / "wgbs_runs.txt")

    genome = resolve_path(project_root, config["reference"]["genome_fasta"])
    gtf = resolve_path(project_root, config["reference"]["gtf"])
    transcripts = resolve_path(project_root, config["reference"]["transcripts_fasta"])
    refs_ready = all(path.exists() for path in (genome, gtf, transcripts))
    docker_ready = shutil.which("docker") is not None
    gpu_enabled = bool(config["parabricks"]["enable_gpu_surrogate_alignment"])

    targets = []
    for run in complete_runs:
        targets.extend(
            [
                f"{outdir}/fastq/trimmed/{run}_1.trim.fastq.gz",
                f"{outdir}/fastq/trimmed/{run}_2.trim.fastq.gz",
                f"{outdir}/qc/fastp/{run}.html",
                f"{outdir}/qc/fastp/{run}.json",
            ]
        )

    if refs_ready and gpu_enabled and docker_ready:
        for run in complete_runs:
            targets.append(f"{outdir}/gpu_surrogate/{run}.fq2bam.bam")
    elif gpu_enabled and complete_runs:
        if not refs_ready:
            warnings.append("Skipping GPU surrogate targets because required reference files are missing.")
        elif not docker_ready:
            warnings.append("Skipping GPU surrogate targets because docker is not installed.")

    if refs_ready and rnaseq_runs:
        for run in complete_runs:
            if run in rnaseq_runs:
                targets.append(f"{outdir}/rnaseq/salmon/{run}/quant.sf")
    elif complete_runs:
        warnings.append("Skipping RNA-seq quantification because references or run metadata are unavailable.")

    if refs_ready and wgbs_runs:
        for run in complete_runs:
            if run in wgbs_runs:
                targets.append(f"{outdir}/wgbs/methyl/{run}/{run}.bismark.cov.gz")
    elif complete_runs:
        warnings.append("Skipping WGBS methylation processing because references or run metadata are unavailable.")

    unclassified = [run for run in complete_runs if run not in rnaseq_runs and run not in wgbs_runs]
    if unclassified:
        warnings.append(f"Skipping assay-specific processing for unclassified runs: {', '.join(unclassified)}")

    return sorted(set(targets)), warnings


def snakemake_cmd(project_root: Path, configfile: Path, cores: int, jobs: int, profile: str, dry_run: bool, targets: list[str]) -> list[str]:
    if shutil.which("snakemake"):
        cmd = ["snakemake"]
    else:
        cmd = ["conda", "run", "-n", "hf_metab_snakemake", "snakemake"]

    cmd.extend(
        [
            "--snakefile",
            "workflow/Snakefile",
            "--configfile",
            str(configfile),
            "--use-conda",
            "--conda-frontend",
            "mamba",
            "--cores",
            str(cores),
            "--jobs",
            str(jobs),
            "--printshellcmds",
            "--rerun-incomplete",
            "--keep-going",
            "--latency-wait",
            "90",
            "--profile",
            profile,
        ]
    )
    if dry_run:
        cmd.append("--dry-run")
    cmd.extend(targets)
    return cmd


def run_once(args: argparse.Namespace, config: dict) -> int:
    raw_dir = resolve_path(args.project_root, config["resources"]["outdir"]) / "fastq" / "raw"
    raw_dir.mkdir(parents=True, exist_ok=True)

    complete_runs, incomplete_runs = discover_complete_runs(raw_dir)
    if incomplete_runs:
        print(f"Incomplete raw pairs still downloading: {', '.join(incomplete_runs)}")
    if not complete_runs:
        print("No complete paired-end raw FASTQ samples are ready yet.")
        return 0

    print(f"Complete paired-end samples ready: {', '.join(complete_runs)}")
    targets, warnings = build_targets(args.project_root, config, complete_runs)
    for warning in warnings:
        print(warning, file=sys.stderr)

    if not targets:
        print("No runnable per-sample targets were identified for the available FASTQs.")
        return 0

    print("Launching incremental per-sample Snakemake targets:")
    for target in targets:
        print(f"  {target}")

    cmd = snakemake_cmd(
        args.project_root,
        args.configfile,
        args.cores,
        args.jobs,
        args.profile,
        args.dry_run,
        targets,
    )
    return subprocess.run(cmd, cwd=args.project_root).returncode


def main() -> int:
    args = parse_args()
    config = load_config(args.configfile)

    if args.watch_interval <= 0:
        return run_once(args, config)

    print(f"Watching for complete FASTQ pairs every {args.watch_interval} seconds.")
    while True:
        exit_code = run_once(args, config)
        if exit_code != 0:
            print(f"Incremental sample run exited with code {exit_code}.", file=sys.stderr)
        time.sleep(args.watch_interval)


if __name__ == "__main__":
    raise SystemExit(main())