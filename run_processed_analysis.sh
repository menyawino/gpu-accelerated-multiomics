#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CORES="${CORES:-4}"
JOBS="${JOBS:-2}"
CONFIGFILE="${PROJECT_ROOT}/config/processed_analysis.yaml"
DRYRUN=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --cores)
      CORES="$2"
      shift 2
      ;;
    --jobs)
      JOBS="$2"
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
    *)
      echo "Unknown argument: $1" >&2
      exit 1
      ;;
  esac
done

CMD=(snakemake -s "${PROJECT_ROOT}/workflow/processed/Snakefile" --configfile "${CONFIGFILE}" --cores "${CORES}" --jobs "${JOBS}" --printshellcmds --rerun-incomplete)
if [[ "${DRYRUN}" -eq 1 ]]; then
  CMD+=(--dry-run)
fi

"${CMD[@]}"
