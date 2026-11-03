#!/bin/bash
# Moves one species/protocol's plots out of a flat benchmark/plot/ folder into
# its own subfolder, e.g. benchmark/plot/ciri_long/.
#
# Every plot file ends in "_<species>" or "_<species>_<strategy suffix>.png",
# where strategy suffix is one of exon_based/strict_bsj/relaxed_bsj/exon/bsj.
# Some plot names themselves contain other species names as a substring
# (for example "classification_circfl_style_ciri_long_exon.png" contains
# "circfl" even for a ciri_long plot), so this strips the known strategy
# suffix first, then checks the species name only at the very end of what
# is left, not as a substring anywhere in the file name.
#
# Usage: sort_plots_by_protocol.sh <plot_dir> <species>
set -euo pipefail

PLOT_DIR="$1"
SPECIES="$2"
SUFFIXES=("_exon_based" "_strict_bsj" "_relaxed_bsj" "_exon" "_bsj")

mkdir -p "$PLOT_DIR/$SPECIES"

shopt -s nullglob
moved=0
for f in "$PLOT_DIR"/*.png; do
    base=$(basename "$f" .png)
    stripped="$base"
    for suf in "${SUFFIXES[@]}"; do
        if [[ "$base" == *"$suf" ]]; then
            stripped="${base%"$suf"}"
            break
        fi
    done
    if [[ "$stripped" == *"_${SPECIES}" ]]; then
        mv "$f" "$PLOT_DIR/$SPECIES/"
        moved=$((moved + 1))
    fi
done
echo "moved $moved plot(s) into $PLOT_DIR/$SPECIES/"
