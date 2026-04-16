#!/usr/bin/env bash
set -euo pipefail

PARABRICKS_IMAGE="${PARABRICKS_IMAGE:-nvcr.io/nvidia/clara/clara-parabricks:4.3.1-1}"
WORKDIR="${PWD}"

if [[ $# -lt 1 ]]; then
  echo "Usage: pbrun.sh <pbrun-subcommand> [args...]" >&2
  exit 1
fi

# Requires prior: docker login nvcr.io
exec docker run --rm --gpus all \
  -u "$(id -u):$(id -g)" \
  -v "${WORKDIR}:${WORKDIR}" \
  -w "${WORKDIR}" \
  "$PARABRICKS_IMAGE" \
  pbrun "$@"
