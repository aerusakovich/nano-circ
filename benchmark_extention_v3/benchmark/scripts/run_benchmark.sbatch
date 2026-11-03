#!/bin/bash
#SBATCH --job-name=nanocirc-benchmark
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --partition=main
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G

# Runs the full benchmark (performance + expression + sensitivity) for one
# species against an nf-core/nanocirc benchmark dataset_root.
#
# Required env vars:
#   DATASET_ROOT   benchmark dataset root (fastq/, ground_truth/, results/, benchmark/)
#   SPECIES        species/dataset key, e.g. mouse or human
#
# Optional: CONFIG (YAML overriding configs/default.yaml), PROJECT_ROOT

set -euo pipefail

: "${DATASET_ROOT:?must set DATASET_ROOT}"
: "${SPECIES:?must set SPECIES}"

PROJECT_ROOT="${PROJECT_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
CONFIG="${CONFIG:-$PROJECT_ROOT/configs/default.yaml}"

python -m nanocirc_benchmark.cli all \
    --dataset-root "$DATASET_ROOT" \
    --species "$SPECIES" \
    --config "$CONFIG"
