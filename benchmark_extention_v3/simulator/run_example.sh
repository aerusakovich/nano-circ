#!/bin/bash
# Quick test: builds a tiny synthetic genome, runs `simulate`, and
# checks that it works. No downloads, no real data needed. Takes seconds.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"

DATA_DIR="example_data"
OUT_DIR="example_results"

echo "=== 1. Building toy genome + GTF (if not already present) ==="
if [ ! -f "$DATA_DIR/toy_genome.fa" ]; then
    python3 scripts/make_toy_data.py --outdir "$DATA_DIR"
else
    echo "Already exists: $DATA_DIR/toy_genome.fa"
fi

echo "=== 2. Running nanocirc-sim simulate ==="
PYTHONPATH=. python3 -m nanocirc_sim.cli simulate \
    --genome "$DATA_DIR/toy_genome.fa" \
    --gtf "$DATA_DIR/toy_annotation.gtf" \
    --total 20 --total-reads 2000 --seed 42 \
    --outdir "$OUT_DIR"

echo "=== Done. Outputs in $OUT_DIR/ ==="
ls "$OUT_DIR"
