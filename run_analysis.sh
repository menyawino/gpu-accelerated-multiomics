#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SNAKE_CORES="${SNAKE_CORES:-8}"
SNAKE_JOBS="${SNAKE_JOBS:-4}"
PROFILE="${PROFILE:-profiles/local}"
ALLOW_NO_GPU="${ALLOW_NO_GPU:-0}"

usage() {
  cat <<'EOF'
Usage: ./run_analysis.sh [--cores N] [--jobs N] [--configfile path] [--dry-run] [--available-only] [--watch-interval N]

This script:
1) Runs setup + dependency checks
2) Launches the Snakemake multi-omics workflow for GSE123976
EOF
}

CONFIGFILE="${PROJECT_ROOT}/config/config.yaml"
DRYRUN=0
AVAILABLE_ONLY=0
WATCH_INTERVAL=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --cores)
      SNAKE_CORES="$2"
      shift 2
      ;;
    --jobs)
      SNAKE_JOBS="$2"
      shift 2
      ;;
    --configfile)
      CONFIGFILE="$2"
      shift 2
      ;;
    --dry-run)
      DRYRUN=1
      shift
      ;;
    --available-only)
      AVAILABLE_ONLY=1
      shift
      ;;
    --watch-interval)
      WATCH_INTERVAL="$2"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Unknown argument: $1" >&2
      usage
      exit 1
      ;;
  esac
done

cd "$PROJECT_ROOT"

if [[ $WATCH_INTERVAL -gt 0 && $AVAILABLE_ONLY -ne 1 ]]; then
  echo "--watch-interval requires --available-only" >&2
  exit 1
fi

SETUP_ARGS=()
if [[ "$ALLOW_NO_GPU" == "1" || $DRYRUN -eq 1 || $AVAILABLE_ONLY -eq 1 ]]; then
  SETUP_ARGS+=(--allow-no-gpu)
fi

bash "$PROJECT_ROOT/setup_tools.sh" "${SETUP_ARGS[@]}"

if [[ $AVAILABLE_ONLY -eq 1 ]]; then
  AVAILABLE_CMD=(
    conda run -n hf_metab python scripts/run_available_samples.py
    --project-root "$PROJECT_ROOT"
    --configfile "$CONFIGFILE"
    --cores "$SNAKE_CORES"
    --jobs "$SNAKE_JOBS"
    --profile "$PROFILE"
  )

  if [[ $DRYRUN -eq 1 ]]; then
    AVAILABLE_CMD+=(--dry-run)
  fi

  if [[ $WATCH_INTERVAL -gt 0 ]]; then
    AVAILABLE_CMD+=(--watch-interval "$WATCH_INTERVAL")
  fi

  "${AVAILABLE_CMD[@]}"
  exit $?
fi

# Prefer global snakemake if available, else use the dedicated conda env.
if command -v snakemake >/dev/null 2>&1; then
  SNAKE_CMD=(snakemake)
else
  SNAKE_CMD=(conda run -n hf_metab_snakemake snakemake)
fi

COMMON_ARGS=(
  --snakefile workflow/Snakefile
  --configfile "$CONFIGFILE"
  --use-conda
  --conda-frontend mamba
  --cores "$SNAKE_CORES"
  --jobs "$SNAKE_JOBS"
  --printshellcmds
  --rerun-incomplete
  --keep-going
  --latency-wait 90
  --profile "$PROFILE"
)

if [[ $DRYRUN -eq 1 ]]; then
  "${SNAKE_CMD[@]}" "${COMMON_ARGS[@]}" --dry-run
else
  "${SNAKE_CMD[@]}" "${COMMON_ARGS[@]}"
fi
