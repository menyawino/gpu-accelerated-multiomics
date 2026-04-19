#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONDA_ENV_NAME="hf_metab"
PARABRICKS_IMAGE="nvcr.io/nvidia/clara/clara-parabricks:4.3.1-1"
CHECK_ONLY=0
ALLOW_NO_GPU=0

log() {
  printf "[%s] %s\n" "$(date +"%Y-%m-%d %H:%M:%S")" "$*"
}

usage() {
  cat <<'EOF'
Usage: ./setup_tools.sh [--check-only] [--allow-no-gpu]

Installs and validates tools needed by the Snakemake multi-omics workflow.

Options:
  --check-only   Only verify required tools and exit
  --allow-no-gpu Continue without NVIDIA validation (useful for local dry-runs)
EOF
}

for arg in "$@"; do
  case "$arg" in
    --check-only)
      CHECK_ONLY=1
      shift
      ;;
    --allow-no-gpu)
      ALLOW_NO_GPU=1
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      printf "Unknown argument: %s\n" "$arg" >&2
      usage
      exit 1
      ;;
  esac
done

require_cmd() {
  local cmd="$1"
  if ! command -v "$cmd" >/dev/null 2>&1; then
    printf "Missing command: %s\n" "$cmd" >&2
    return 1
  fi
}

install_miniforge_if_missing() {
  if command -v conda >/dev/null 2>&1; then
    return 0
  fi

  local os arch installer url target
  os="$(uname -s)"
  arch="$(uname -m)"

  case "$os" in
    Linux)
      case "$arch" in
        x86_64) installer="Miniforge3-Linux-x86_64.sh" ;;
        aarch64|arm64) installer="Miniforge3-Linux-aarch64.sh" ;;
        *)
          printf "Unsupported Linux architecture: %s\n" "$arch" >&2
          return 1
          ;;
      esac
      ;;
    Darwin)
      case "$arch" in
        x86_64) installer="Miniforge3-MacOSX-x86_64.sh" ;;
        arm64) installer="Miniforge3-MacOSX-arm64.sh" ;;
        *)
          printf "Unsupported macOS architecture: %s\n" "$arch" >&2
          return 1
          ;;
      esac
      ;;
    *)
      printf "Unsupported OS for auto-installing conda: %s\n" "$os" >&2
      return 1
      ;;
  esac

  url="https://github.com/conda-forge/miniforge/releases/latest/download/${installer}"
  target="${HOME}/miniforge3"

  log "Installing Miniforge to ${target}"
  curl -L "$url" -o /tmp/miniforge.sh
  bash /tmp/miniforge.sh -b -p "$target"
  rm -f /tmp/miniforge.sh

  # shellcheck disable=SC1091
  source "${target}/etc/profile.d/conda.sh"
  conda init bash >/dev/null 2>&1 || true
}

ensure_conda_shell() {
  if command -v conda >/dev/null 2>&1; then
    local conda_base
    conda_base="$(conda info --base)"
    # shellcheck disable=SC1091
    source "${conda_base}/etc/profile.d/conda.sh"
  fi
}

ensure_mamba() {
  if command -v mamba >/dev/null 2>&1; then
    return 0
  fi

  log "Installing mamba into the base conda environment"
  conda install -n base -c conda-forge -y mamba

  if ! command -v mamba >/dev/null 2>&1; then
    hash -r
  fi

  command -v mamba >/dev/null 2>&1
}

run_conda_env_update() {
  local log_file
  log_file="$(mktemp)"

  if mamba env update -n "$CONDA_ENV_NAME" -f "${PROJECT_ROOT}/envs/workflow.yaml" --prune \
    2>&1 | tee "$log_file"; then
    rm -f "$log_file"
    return 0
  fi

  if grep -Fq "Cannot link a source that does not exist" "$log_file"; then
    log "Package cache appears stale; cleaning cache and retrying ${CONDA_ENV_NAME} once with mamba"
    conda clean --packages --tarballs --yes
    conda env remove -n "$CONDA_ENV_NAME" --yes >/dev/null 2>&1 || true
    mamba env update -n "$CONDA_ENV_NAME" -f "${PROJECT_ROOT}/envs/workflow.yaml" --prune
    rm -f "$log_file"
    return 0
  fi

  rm -f "$log_file"
  return 1
}

install_workflow_dependencies() {
  ensure_conda_shell
  if ! command -v conda >/dev/null 2>&1; then
    printf "Conda is still unavailable after attempted install.\n" >&2
    return 1
  fi

  ensure_mamba

  log "Creating/updating workflow env with mamba: ${CONDA_ENV_NAME}"
  conda config --set channel_priority strict
  run_conda_env_update

  log "Installing Snakemake launcher env with mamba (if needed)"
  if ! command -v snakemake >/dev/null 2>&1; then
    mamba create -y -n hf_metab_snakemake -c conda-forge -c bioconda snakemake=8.30.0 mamba
  fi
}

write_parabricks_wrapper() {
  local wrapper
  wrapper="${PROJECT_ROOT}/scripts/pbrun.sh"

  cat > "$wrapper" <<'EOF'
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
EOF
  chmod +x "$wrapper"
}

validate_gpu_stack() {
  local missing=0
  require_cmd docker || missing=1

  if command -v nvidia-smi >/dev/null 2>&1; then
    nvidia-smi >/dev/null 2>&1 || {
      printf "nvidia-smi exists but failed to query GPU.\n" >&2
      missing=1
    }
  else
    printf "nvidia-smi not found. GPU steps require an NVIDIA Linux host with drivers.\n" >&2
    missing=1
  fi

  if [[ $missing -eq 1 ]]; then
    return 1
  fi

  log "GPU stack check passed"
}

validate_workflow_tools() {
  local missing=0
  local cmds=(python3 samtools fastp prefetch fasterq-dump STAR salmon bismark bismark_methylation_extractor bedtools)

  for cmd in "${cmds[@]}"; do
    if ! command -v "$cmd" >/dev/null 2>&1; then
      printf "Missing command in PATH: %s\n" "$cmd" >&2
      missing=1
    fi
  done

  if ! command -v snakemake >/dev/null 2>&1; then
    log "snakemake not in current PATH. Will still run via 'conda run -n hf_metab_snakemake snakemake' in analysis script."
  fi

  if [[ $missing -eq 1 ]]; then
    return 1
  fi

  log "Workflow tools check passed"
}

main() {
  log "Project root: ${PROJECT_ROOT}"

  if [[ $CHECK_ONLY -eq 0 ]]; then
    install_miniforge_if_missing
    install_workflow_dependencies
    write_parabricks_wrapper
  fi

  ensure_conda_shell

  if [[ -d "$(conda info --base 2>/dev/null)/envs/${CONDA_ENV_NAME}" ]]; then
    # shellcheck disable=SC1091
    conda activate "$CONDA_ENV_NAME"
  fi

  validate_workflow_tools

  if ! validate_gpu_stack; then
    if [[ $ALLOW_NO_GPU -eq 1 ]]; then
      log "GPU validation failed but --allow-no-gpu was set; continuing."
      exit 0
    fi
    cat <<'EOF'
GPU validation failed.
If you are currently on macOS for development, run analysis on a Linux NVIDIA host (e.g., T4 VM) with:
1) NVIDIA driver + container toolkit
2) Docker daemon configured for --gpus all
3) Access to nvcr.io for Parabricks image pulls
EOF
    exit 1
  fi

  log "Setup and validation completed successfully"
}

main "$@"
