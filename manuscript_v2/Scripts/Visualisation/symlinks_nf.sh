#!/usr/bin/env bash
set -u

WORK_DIR="/projects/dog/anastasia/PCI_benchmark/run1_mouse_in_silico_200000/work" # Adjust for each run
OUTPUT_DIR="/projects/dog/anastasia/PCI_benchmark/run1_mouse_in_silico_200000/visualisation" # Adjust for each run
DRY_RUN=false
MAX_LEVELS_UP=8
MAX_RUN=10

TARGET_FILES=(
  "CIRI-long.cand_circ.bed12"
  "comprehensive_circrna.bed12"
  "ground_truth_circular.bed12"
  "ground_truth.bed"
  "isocirc.bed"
  "CIRI-long.expression"
  "combined_reads.circRNA_candidates.annotated.txt"
  "time_vs_memory_performance.png"
  "throughput_vs_memory_performance.png"
  "isocirc.out"
)

# Files that are "global" (not per-run, no run_X association expected)
GLOBAL_FILES=(
  "ground_truth.bed"
  "time_vs_memory_performance.png"
  "throughput_vs_memory_performance.png"
)

log() { printf "%s\n" "$*" >&2; }

TEMP="$(mktemp -d)"
trap 'rm -rf "$TEMP"' EXIT

ALL="$TEMP/all_files.tsv"
: > "$ALL"

abspath() { readlink -f "$1" 2>/dev/null || true; }

is_in_workdir() {
  local p="$1"
  [[ -n "$p" && "$p" == "$WORK_DIR"* ]]
}

find_run_from_dir_tree() {
  local start_dir="$1"
  local d="$start_dir"
  local level i r tgt

  for ((level=0; level<MAX_LEVELS_UP; level++)); do
    for ((i=1; i<=MAX_RUN; i++)); do
      r="$d/run_$i"
      if [[ -L "$r" ]]; then
        tgt="$(abspath "$r")"
        if [[ -d "$tgt" ]] && is_in_workdir "$tgt"; then
          printf "%s|%s\n" "$i" "$tgt"
          return 0
        fi
      fi
    done
    [[ "$d" == "$WORK_DIR" || "$d" == "/" ]] && break
    d="$(dirname "$d")"
  done
  return 1
}

determine_run_for_file() {
  local filename="$1"
  local filepath="$2"

  local resolved="$filepath"
  if [[ -L "$filepath" ]]; then
    local rf
    rf="$(abspath "$filepath")"
    [[ -n "$rf" ]] && resolved="$rf"
  fi
  local filedir
  filedir="$(dirname "$resolved")"

  # Tool-mapped files
  if [[ "$filename" == "CIRI-long.cand_circ.bed12" ]]; then
    if [[ -L "$filedir/cirilong_output" ]]; then
      local t tdir rr
      t="$(abspath "$filedir/cirilong_output")"
      [[ -n "$t" ]] || return 1
      tdir="$(dirname "$t")"
      rr="$(find_run_from_dir_tree "$tdir" 2>/dev/null)" || return 1
      printf "%s|%s|%s\n" "${rr%%|*}" "${rr#*|}" "$resolved"
      return 0
    fi
  fi

  if [[ "$filename" == "comprehensive_circrna.bed12" ]]; then
    if [[ -L "$filedir/circnick_output" ]]; then
      local t tdir rr
      t="$(abspath "$filedir/circnick_output")"
      [[ -n "$t" ]] || return 1
      tdir="$(dirname "$t")"
      rr="$(find_run_from_dir_tree "$tdir" 2>/dev/null)" || return 1
      printf "%s|%s|%s\n" "${rr%%|*}" "${rr#*|}" "$resolved"
      return 0
    fi
  fi

  # ground_truth_circular
  if [[ "$filename" == "ground_truth_circular.bed12" ]]; then
    local rr
    rr="$(find_run_from_dir_tree "$filedir" 2>/dev/null)" || true
    if [[ -n "${rr:-}" ]]; then
      printf "%s|%s|%s\n" "${rr%%|*}" "${rr#*|}" "$resolved"
      return 0
    fi
  fi

  # Generic
  local rr
  rr="$(find_run_from_dir_tree "$filedir" 2>/dev/null)" || return 1
  printf "%s|%s|%s\n" "${rr%%|*}" "${rr#*|}" "$resolved"
  return 0
}

echo "========================================================================"
echo "NEXTFLOW ORGANIZER"
echo "========================================================================"
echo ""

log "Searching for files under: $WORK_DIR"
log ""

for filename in "${TARGET_FILES[@]}"; do
  log "Processing: $filename"

  while IFS= read -r -d '' filepath; do
    # Skip globals here (handled later)
    for g in "${GLOBAL_FILES[@]}"; do
      [[ "$filename" == "$g" ]] && continue 2
    done

    out="$(determine_run_for_file "$filename" "$filepath" 2>/dev/null)" || continue
    run_num="${out%%|*}"
    rest="${out#*|}"
    run_target="${rest%%|*}"
    resolved_file="${rest#*|}"

    ts="$(stat -c %Y "$run_target" 2>/dev/null || stat -f %m "$run_target" 2>/dev/null || echo 0)"
    printf "%s\t%s\t%s\t%s\t%s\t%s\n" \
      "$ts" "$run_num" "$run_target" "$filename" "$filepath" "$resolved_file" >> "$ALL"
  done < <(find "$WORK_DIR" -name "$filename" -print0 2>/dev/null)

  count="$(awk -F'\t' -v f="$filename" '$4==f{c++} END{print c+0}' "$ALL" 2>/dev/null || echo 0)"
  log "  → Found $count with TRUSTED run associations"
done

if [[ ! -s "$ALL" ]]; then
  echo ""
  echo "ERROR: No files found with trusted run associations."
  exit 1
fi

echo ""
echo "========================================================================"
echo "GROUPING BY (run_num, run_target)"
echo "========================================================================"
echo ""

GROUPS="$TEMP/groups.tsv"
# Build groups without fragile pipelines
awk -F'\t' '{print $1"\t"$2"\t"$3}' "$ALL" | sort -u | sort -k1,1n -k2,2n > "$GROUPS"

if [[ ! -s "$GROUPS" ]]; then
  echo "ERROR: GROUPS is empty. Showing first 20 lines of ALL:"
  head -n 20 "$ALL"
  exit 1
fi

# Print batches
batch=0
while IFS=$'\t' read -r ts run_num run_target; do
  ((batch++))
  nfiles="$(awk -F'\t' -v rn="$run_num" -v rt="$run_target" '$2==rn && $3==rt{c++} END{print c+0}' "$ALL")"

  echo "Batch $batch: run_${run_num}"
  echo "  Run link target: $run_target"
  if command -v date >/dev/null 2>&1; then
    if date -d "@$ts" "+%Y-%m-%d %H:%M" >/dev/null 2>&1; then
      echo "  Time: $(date -d "@$ts" "+%Y-%m-%d %H:%M")"
    else
      echo "  Time: $ts"
    fi
  fi
  echo "  Files: $nfiles"
  awk -F'\t' -v rn="$run_num" -v rt="$run_target" '$2==rn && $3==rt{print $4}' "$ALL" | sort -u | sed 's/^/    /'
  echo ""
done < "$GROUPS"

echo "========================================================================"
echo "OUTPUT STRUCTURE"
echo "========================================================================"
echo ""

batch=0
while IFS=$'\t' read -r ts run_num run_target; do
  ((batch++))
  nfiles="$(awk -F'\t' -v rn="$run_num" -v rt="$run_target" '$2==rn && $3==rt{c++} END{print c+0}' "$ALL")"
  echo "not_nf/batch_${batch}/run_${run_num}/ ($nfiles files + globals)"
  awk -F'\t' -v rn="$run_num" -v rt="$run_target" '$2==rn && $3==rt{print $4}' "$ALL" | sort -u | sed 's/^/  /'
  printf "  %s\n" "${GLOBAL_FILES[@]}"
  echo ""
done < "$GROUPS"

if [[ "$DRY_RUN" == "true" ]]; then
  echo "========================================================================"
  echo "DRY RUN - Set DRY_RUN=false to copy files"
  exit 0
fi

echo "========================================================================"
read -r -p "Copy files? (yes/no): " ans
[[ "$ans" != "yes" && "$ans" != "y" ]] && exit 0

echo ""
echo "Copying..."
mkdir -p "$OUTPUT_DIR"
MAPPING="$OUTPUT_DIR/batch_mapping.txt"
: > "$MAPPING"

# Locate global files once (anywhere under WORK_DIR, take first hit)
declare -A GLOBAL_PATH
for g in "${GLOBAL_FILES[@]}"; do
  gp="$(find "$WORK_DIR" -name "$g" -print -quit 2>/dev/null || true)"
  GLOBAL_PATH["$g"]="$gp"
done

batch=0
while IFS=$'\t' read -r ts run_num run_target; do
  ((batch++))
  dest="$OUTPUT_DIR/batch_${batch}/run_${run_num}"
  mkdir -p "$dest"

  echo ""
  echo "batch_${batch}/run_${run_num}/"

  {
    echo "batch_${batch}/run_${run_num}/"
    echo "  Run link target: $run_target"
  } >> "$MAPPING"

  # Copy run-associated files
  awk -F'\t' -v rn="$run_num" -v rt="$run_target" '$2==rn && $3==rt' "$ALL" \
  | while IFS=$'\t' read -r _ts _rn _rt fname orig resolved; do
      src="$resolved"
      [[ -f "$src" ]] || src="$orig"
      if [[ -f "$src" ]]; then
        cp -p "$src" "$dest/$fname"
        echo "  ✓ $fname"
        {
          echo "    $fname"
          echo "      Source: $src"
        } >> "$MAPPING"
      else
        echo "  ✗ $fname"
        {
          echo "    $fname"
          echo "      Source: MISSING (orig=$orig resolved=$resolved)"
        } >> "$MAPPING"
      fi
    done

  # Copy globals into every batch
  for g in "${GLOBAL_FILES[@]}"; do
    gp="${GLOBAL_PATH[$g]}"
    if [[ -n "$gp" && -f "$gp" ]]; then
      cp -p "$gp" "$dest/$g"
      echo "  ✓ $g (global)"
    else
      echo "  ✗ $g (global not found)"
    fi
  done

  echo "" >> "$MAPPING"
done < "$GROUPS"

echo ""
echo "DONE! Output: $OUTPUT_DIR"
echo "Mapping: $MAPPING"
