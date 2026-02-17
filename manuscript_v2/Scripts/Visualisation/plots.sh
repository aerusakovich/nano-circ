#!/usr/bin/env bash

# ------------------------------------------------------------------------------
# ENVIRONMENT SETUP
# ------------------------------------------------------------------------------
source /local/env/envpython-3.9.5.sh
conda activate bed12

if ! command -v Rscript &>/dev/null; then
  if [ -f "/local/env/envR-4.0.3.sh" ]; then
    source /local/env/envR-4.0.3.sh
  fi
fi
set -euo pipefail
# ------------------------------------------------------------------------------
# CONFIGURATION
# ------------------------------------------------------------------------------
BASE_DIR="/projects/dog/anastasia/PCI_benchmark/run1_mouse_in_silico_200000/visualisation" # Adjust as needed
OUTPUT_DIR="${BASE_DIR}/clean_plots" # Adjust as needed
GTF_FILE="/home/genouest/cnrs_umr6290/aerusakovich/benchmark/data/ciri-long/files/m38/gencode.vM10.annotation.gtf"

# Analysis parameters
FRACTIONS=(0.25 0.5 0.75 0.95) # Adjust as needed
SPLIT_MODES=(true false)

# File names per tool
GT_BED_NAME="ground_truth_circular.bed12"
CIRI_BED_NAME="CIRI-long.cand_circ.bed12"
ISO_BED_NAME="isocirc.bed"
NICK_BED_NAME="comprehensive_circrna.bed12"

# Expression file names
CIRI_EXPR_NAME="CIRI-long.expression"
ISO_EXPR_NAME="isocirc.out"
NICK_EXPR_NAME="combined_reads.circRNA_candidates.annotated.txt"

# Create output structure
mkdir -p "${OUTPUT_DIR}"/{data,plots}
mkdir -p "${OUTPUT_DIR}/data"/{sanitized,length,classification,performance,performance_filtered,combos,combos_filtered,expression,upset}
mkdir -p "${OUTPUT_DIR}/plots"/{length,classification,performance,performance_filtered,combos,combos_filtered,expression,upset}

echo "=============================================================================="
echo "CircRNA Analysis Pipeline - Clean Version"
echo "=============================================================================="
echo "Base directory:   ${BASE_DIR}"
echo "Output directory: ${OUTPUT_DIR}"
echo "GTF file:         ${GTF_FILE}"
echo ""

# ------------------------------------------------------------------------------
# UTILITY FUNCTIONS
# ------------------------------------------------------------------------------
run_dirs() {
  find "${BASE_DIR}" -maxdepth 1 -type d -name "run_*" | sort
}

log_step() {
  echo ""
  echo "=============================================================================="
  echo "$1"
  echo "=============================================================================="
}

log_substep() {
  echo "  → $1"
}

needs_run() {
  local output_file="$1"
  if [ -f "$output_file" ]; then
    echo "  ✓ Already exists: $(basename "$output_file")"
    return 1
  fi
  return 0
}

################################################################################
# STEP 1: DATA SANITIZATION 
################################################################################
log_step "STEP 1: Data Sanitization"

SANITIZED_DIR="${OUTPUT_DIR}/data/sanitized"
mkdir -p "${SANITIZED_DIR}"

SANITIZE_MARKER="${SANITIZED_DIR}/.sanitization_complete"

if [ -f "$SANITIZE_MARKER" ]; then
  echo "  ✓ Sanitization already complete"
else
  log_substep "Sanitizing BED12 files (removing malformed records)"
  
  sanitize_bed12() {
    local input_bed="$1"
    local output_bed="$2"
    local label="$3"
    
    [ -f "$input_bed" ] || { echo "    $label: FILE NOT FOUND"; return 1; }
    
    local raw_count=$(wc -l < "$input_bed")
    
    # More permissive sanitization - only remove severely malformed records
    awk -F'\t' 'BEGIN{OFS="\t"}
      /^#/ || /^track/ || /^browser/ { next }
      {
        # Check minimum requirements
        if (NF < 12) next;
        
        s=$2; e=$3;
        # Allow equal start/end (single position), just not invalid ranges
        if (s=="" || e=="" || s<0 || e<0 || s>e) next;
        
        # Only check that blockSizes and blockStarts exist, dont validate values strictly
        bs=$11; bo=$12;
        if (bs=="" || bo=="") next;
        
        # Check for obviously bad values (negative)
        nbs=split(bs,a,","); 
        bad=0;
        for(i=1; i<=nbs; i++) {
          if (a[i]!="" && a[i]<0) { bad=1; break }
        }
        if(bad) next;
        
        nbo=split(bo,b,",");
        for(i=1; i<=nbo; i++) {
          if (b[i]!="" && b[i]<0) { bad=1; break }
        }
        if(bad) next;
        
        # Passed all checks
        print $0;
      }' "$input_bed" > "$output_bed"
    
    local clean_count=$(wc -l < "$output_bed")
    local dropped=$((raw_count - clean_count))
    
    if [ "$clean_count" -eq 0 ]; then
      echo "    ERROR: $label sanitization resulted in EMPTY file! Using original instead."
      cp "$input_bed" "$output_bed"
      clean_count=$raw_count
      dropped=0
    fi
    
    echo "    $label: raw=$raw_count, clean=$clean_count, dropped=$dropped"
    
    # Verify output is readable by bedtools
    if ! grep -v '^#' "$output_bed" | head -1 | awk -F'\t' '{if(NF<3) exit 1}'; then
      echo "    ERROR: $label output appears malformed! Using original."
      cp "$input_bed" "$output_bed"
    fi
    
    return 0
  }
  
  for rd in $(run_dirs); do
    run_name="$(basename "$rd")"
    log_substep "Processing ${run_name}"
    
    run_sanitized_dir="${SANITIZED_DIR}/${run_name}"
    mkdir -p "$run_sanitized_dir"
    
    # Check if source files exist
    echo "    Checking source files..."
    [ -f "${rd}/${GT_BED_NAME}" ] && echo "      ✓ GT found" || echo "      ✗ GT MISSING"
    [ -f "${rd}/${CIRI_BED_NAME}" ] && echo "      ✓ CIRI found" || echo "      ✗ CIRI missing"
    [ -f "${rd}/${ISO_BED_NAME}" ] && echo "      ✓ IsoCirc found" || echo "      ✗ IsoCirc missing"
    [ -f "${rd}/${NICK_BED_NAME}" ] && echo "      ✓ CircNick found" || echo "      ✗ CircNick missing"
    
    sanitize_bed12 "${rd}/${GT_BED_NAME}" "${run_sanitized_dir}/${GT_BED_NAME}" "GT"
    sanitize_bed12 "${rd}/${CIRI_BED_NAME}" "${run_sanitized_dir}/${CIRI_BED_NAME}" "CIRI"
    sanitize_bed12 "${rd}/${ISO_BED_NAME}" "${run_sanitized_dir}/${ISO_BED_NAME}" "IsoCirc"
    sanitize_bed12 "${rd}/${NICK_BED_NAME}" "${run_sanitized_dir}/${NICK_BED_NAME}" "CircNick"
    
    # Verify sanitized files
    echo "    Verifying sanitized files..."
    for f in "${GT_BED_NAME}" "${CIRI_BED_NAME}" "${ISO_BED_NAME}" "${NICK_BED_NAME}"; do
      if [ -f "${run_sanitized_dir}/${f}" ]; then
        count=$(wc -l < "${run_sanitized_dir}/${f}")
        echo "      ${f}: ${count} lines"
      fi
    done
    
    # Copy expression files (no sanitization needed)
    [ -f "${rd}/${CIRI_EXPR_NAME}" ] && cp "${rd}/${CIRI_EXPR_NAME}" "${run_sanitized_dir}/"
    [ -f "${rd}/${ISO_EXPR_NAME}" ] && cp "${rd}/${ISO_EXPR_NAME}" "${run_sanitized_dir}/"
    [ -f "${rd}/${NICK_EXPR_NAME}" ] && cp "${rd}/${NICK_EXPR_NAME}" "${run_sanitized_dir}/"
  done
  
  touch "$SANITIZE_MARKER"
  echo "  ✓ Sanitization complete"
fi

# Prepare GTF features (run once)
log_step "STEP 1B: GTF Feature Extraction"

GTF_FEATURES_DIR="${OUTPUT_DIR}/data/gtf_features"
mkdir -p "$GTF_FEATURES_DIR"

if needs_run "${GTF_FEATURES_DIR}/genes.bed"; then
  log_substep "Extracting gene features"
  
  awk '$3 == "gene"' "$GTF_FILE" | \
    awk 'BEGIN{OFS="\t"} {print $1, $4-1, $5, ".", 0, $7}' > "${GTF_FEATURES_DIR}/genes.bed"
  
  echo "  ✓ Gene features extracted: $(wc -l < "${GTF_FEATURES_DIR}/genes.bed") genes"
fi

if needs_run "${GTF_FEATURES_DIR}/exons.bed"; then
  log_substep "Extracting exon features"
  
  awk '$3 == "exon"' "$GTF_FILE" | \
    awk 'BEGIN{OFS="\t"} {print $1, $4-1, $5, ".", 0, $7}' > "${GTF_FEATURES_DIR}/exons.bed"
  
  echo "  ✓ Exon features extracted: $(wc -l < "${GTF_FEATURES_DIR}/exons.bed") exons"
fi

################################################################################
# STEP 2: LENGTH ANALYSIS
################################################################################
log_step "STEP 2: Length Analysis"

LENGTH_DATA_DIR="${OUTPUT_DIR}/data/length"
LENGTH_PLOT_DIR="${OUTPUT_DIR}/plots/length"
mkdir -p "$LENGTH_DATA_DIR" "$LENGTH_PLOT_DIR"

LENGTH_CSV="${LENGTH_DATA_DIR}/length_data.csv"

if needs_run "$LENGTH_CSV"; then
  log_substep "Collecting length data from sanitized BED12 files"
  
  echo "Tool,Run,Length" > "$LENGTH_CSV"
  
  calc_mature_len_from_bed12() {
    local tool_id="$1"
    local run_name="$2"
    local bed12="$3"
    local output_csv="$4"
    
    [ -f "$bed12" ] || return 0
    
    awk -v tool="$tool_id" -v run="$run_name" '
      BEGIN{FS=OFS="\t"}
      !/^#/ && !/^track/ && !/^browser/ {
        n=split($11, sizes, ",");
        len=0;
        for(i=1;i<=n;i++) if(sizes[i]!="") len+=sizes[i];
        print tool "," run "," len
      }' "$bed12" >> "$output_csv"
  }
  
  for rd in $(run_dirs); do
    run_name="$(basename "$rd")"
    sanitized_run="${SANITIZED_DIR}/${run_name}"
    
    calc_mature_len_from_bed12 "ground_truth" "$run_name" "${sanitized_run}/${GT_BED_NAME}" "$LENGTH_CSV"
    calc_mature_len_from_bed12 "ciri-long" "$run_name" "${sanitized_run}/${CIRI_BED_NAME}" "$LENGTH_CSV"
    calc_mature_len_from_bed12 "isocirc" "$run_name" "${sanitized_run}/${ISO_BED_NAME}" "$LENGTH_CSV"
    calc_mature_len_from_bed12 "circnick" "$run_name" "${sanitized_run}/${NICK_BED_NAME}" "$LENGTH_CSV"
  done
  
  echo "  ✓ Collected $(( $(wc -l < "$LENGTH_CSV") - 1 )) length measurements"
fi

# Generate length barplot only
if needs_run "${LENGTH_PLOT_DIR}/length_barplot.png"; then
  log_substep "Generating length barplot"
  
  cd "$LENGTH_PLOT_DIR"
  cp "$LENGTH_CSV" ./
  
  Rscript - <<'RSCRIPT'
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(scales)
})

data <- read.csv("length_data.csv", stringsAsFactors=FALSE)

tool_order  <- c("ground_truth","ciri-long","isocirc","circnick")
tool_labels <- c("Ground Truth","CIRI-long","IsoCirc","CircNick-lrs")
tool_colors <- c("#808080","#0072B2","#009E73","#D55E00")

data$Tool <- factor(data$Tool, levels=tool_order)

stats <- data %>%
  group_by(Tool) %>%
  summarize(Count=n(), Mean=mean(Length, na.rm=TRUE),
            SE=sd(Length, na.rm=TRUE)/sqrt(n()), .groups="drop")
stats$Tool <- factor(stats$Tool, levels=tool_order)

png("length_barplot.png", width=1800, height=1400, res=150)
p <- ggplot(stats, aes(x=Tool, y=Mean, fill=Tool)) +
  geom_bar(stat="identity", color="black", alpha=0.9, show.legend=FALSE) +
  geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE), width=0.3, linewidth=0.8) +
  scale_fill_manual(values=tool_colors) +
  scale_x_discrete(labels=tool_labels) +
  scale_y_continuous(labels=scales::comma, expand=expansion(mult=c(0,0.15))) +
  labs(title="Mean CircRNA Mature Length by Tool",
       subtitle="Error bars show standard error (SE)", x="", y="Mean Mature Length (bp)") +
  theme_minimal(base_size=20) +
  theme(plot.title=element_text(size=24, face="bold", hjust=0.5),
        axis.text.x=element_text(size=20, angle=15, hjust=1)) +
  geom_text(aes(label=paste0(round(Mean,0),"±",round(SE,1))), vjust=-1.5, size=7, fontface="bold") +
  geom_text(aes(y=Mean/2, label=paste0("n=",scales::comma(Count))), size=7, fontface="bold", color="white")
print(p); dev.off()

cat("✓ Length barplot generated\n")
RSCRIPT
  
  rm length_data.csv
  cd - > /dev/null
  echo "  ✓ Length plot complete"
fi

################################################################################
# STEP 3: CLASSIFICATION ANALYSIS
################################################################################
log_step "STEP 3: Classification Analysis"

CLASS_DATA_DIR="${OUTPUT_DIR}/data/classification"
CLASS_PLOT_DIR="${OUTPUT_DIR}/plots/classification"
mkdir -p "$CLASS_DATA_DIR" "$CLASS_PLOT_DIR"

# Run classification for each combination
CLASSIFICATION_MARKER="${CLASS_DATA_DIR}/.classification_complete"

if [ ! -f "$CLASSIFICATION_MARKER" ]; then
  log_substep "Running classification analysis"
  
  for rd in $(run_dirs); do
    run_name="$(basename "$rd")"
    sanitized_run="${SANITIZED_DIR}/${run_name}"
    
    for frac in "${FRACTIONS[@]}"; do
      for split in "${SPLIT_MODES[@]}"; do
        mode_dir="nosplit"
        [ "$split" = "true" ] && mode_dir="split"
        
        out_dir="${CLASS_DATA_DIR}/${run_name}/f${frac}/${mode_dir}"
        result_file="${out_dir}/classification_results.csv"
        
        if [ -f "$result_file" ]; then
          continue
        fi
        
        mkdir -p "$out_dir/temp"
        
        # Run Python classification
        export CLASS_GTF_FEATURES="${GTF_FEATURES_DIR}"
        export CLASS_GT_BED="${sanitized_run}/${GT_BED_NAME}"
        export CLASS_CIRI_BED="${sanitized_run}/${CIRI_BED_NAME}"
        export CLASS_ISO_BED="${sanitized_run}/${ISO_BED_NAME}"
        export CLASS_NICK_BED="${sanitized_run}/${NICK_BED_NAME}"
        export CLASS_RUN_NAME="$run_name"
        export CLASS_FRACTION="$frac"
        export CLASS_SPLIT="$split"
        export CLASS_OUT_DIR="$out_dir"
        
        python3 - <<'PYCLASSIFY'
import os, subprocess
from collections import defaultdict

GTF_FEATURES = os.environ["CLASS_GTF_FEATURES"]
GT_BED = os.environ["CLASS_GT_BED"]
CIRI_BED = os.environ.get("CLASS_CIRI_BED", "")
ISO_BED = os.environ.get("CLASS_ISO_BED", "")
NICK_BED = os.environ.get("CLASS_NICK_BED", "")
RUN_NAME = os.environ["CLASS_RUN_NAME"]
FRACTION = os.environ["CLASS_FRACTION"]
SPLIT = os.environ["CLASS_SPLIT"]
OUT_DIR = os.environ["CLASS_OUT_DIR"]
TEMP_DIR = os.path.join(OUT_DIR, "temp")

gene_bed = os.path.join(GTF_FEATURES, "genes.bed")
exon_bed = os.path.join(GTF_FEATURES, "exons.bed")

def run_cmd(cmd):
    r = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    return r.stdout

def read_coords(bed_file):
    coords = set()
    if not os.path.exists(bed_file): return coords
    with open(bed_file) as f:
        for line in f:
            if line.startswith('#') or line.startswith('track'): continue
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                coords.add((parts[0], int(parts[1]), int(parts[2])))
    return coords

def classify_coords(coords_set, label):
    if not coords_set: return defaultdict(int)
    
    os.makedirs(TEMP_DIR, exist_ok=True)
    tmp_bed = os.path.join(TEMP_DIR, f"{label}.bed")
    
    with open(tmp_bed, 'w') as f:
        for c, s, e in sorted(coords_set):
            f.write(f"{c}\t{s}\t{e}\t.\t0\t+\n")
    
    def extract(cmd):
        out = run_cmd(cmd).strip().split('\n')
        coords = set()
        for line in out:
            if not line.strip(): continue
            p = line.split('\t')
            if len(p) >= 3: coords.add((p[0], int(p[1]), int(p[2])))
        return coords
    
    gene_s = extract(f"bedtools intersect -a {tmp_bed} -b {gene_bed} -s -wa -wb")
    exon_s = extract(f"bedtools intersect -a {tmp_bed} -b {exon_bed} -s -wa -wb")
    full_exon_s = extract(f"bedtools intersect -a {tmp_bed} -b {exon_bed} -s -f 1.0 -wa")
    
    gene_u = extract(f"bedtools intersect -a {tmp_bed} -b {gene_bed} -wa -wb")
    exon_u = extract(f"bedtools intersect -a {tmp_bed} -b {exon_bed} -wa -wb")
    full_exon_u = extract(f"bedtools intersect -a {tmp_bed} -b {exon_bed} -f 1.0 -wa")
    
    counts = defaultdict(int)
    for coord in coords_set:
        in_gene = (coord in gene_s) or (coord in gene_u)
        in_full_exon = (coord in full_exon_s) or (coord in full_exon_u)
        in_exon = (coord in exon_s) or (coord in exon_u)
        
        if in_gene:
            if in_full_exon: typ = "eciRNA"
            elif in_exon: typ = "EIciRNA"
            else: typ = "ciRNA"
        else:
            typ = "intergenic"
        counts[typ] += 1
    
    return counts

# Classify ground truth
gt_coords = read_coords(GT_BED)
gt_counts = classify_coords(gt_coords, "gt")

results = []
results.append({
    'Tool': 'ground_truth',
    'Run': RUN_NAME,
    'Fraction': FRACTION,
    'SplitMode': 'split' if SPLIT == 'true' else 'nosplit',
    'eciRNA': gt_counts['eciRNA'],
    'EIciRNA': gt_counts['EIciRNA'],
    'ciRNA': gt_counts['ciRNA'],
    'intergenic': gt_counts['intergenic']
})

# Classify tool predictions (TP only)
tools = []
if os.path.exists(CIRI_BED): tools.append(('ciri-long', CIRI_BED))
if os.path.exists(ISO_BED): tools.append(('isocirc', ISO_BED))
if os.path.exists(NICK_BED): tools.append(('circnick', NICK_BED))

for tool_name, tool_bed in tools:
    # Get TP coords (intersect with GT)
    split_opt = "-split" if SPLIT == "true" else ""
    tp_file = os.path.join(TEMP_DIR, f"{tool_name}_tp.bed")
    
    cmd = f"bedtools intersect -a {GT_BED} -b {tool_bed} -f {FRACTION} -r -wa -u {split_opt} > {tp_file}"
    os.system(cmd)
    
    tp_coords = read_coords(tp_file)
    tp_counts = classify_coords(tp_coords, f"{tool_name}_tp")
    
    results.append({
        'Tool': tool_name,
        'Run': RUN_NAME,
        'Fraction': FRACTION,
        'SplitMode': 'split' if SPLIT == 'true' else 'nosplit',
        'eciRNA': tp_counts['eciRNA'],
        'EIciRNA': tp_counts['EIciRNA'],
        'ciRNA': tp_counts['ciRNA'],
        'intergenic': tp_counts['intergenic']
    })

# Save results
import pandas as pd
df = pd.DataFrame(results)
df.to_csv(os.path.join(OUT_DIR, "classification_results.csv"), index=False)

print(f"✓ Classification: {RUN_NAME} f={FRACTION} {'split' if SPLIT=='true' else 'nosplit'}")
PYCLASSIFY
        
      done
    done
  done
  
  touch "$CLASSIFICATION_MARKER"
  echo "  ✓ Classification complete"
fi

# Aggregate classification results
ALL_CLASS_CSV="${CLASS_DATA_DIR}/all_classification_results.csv"

if needs_run "$ALL_CLASS_CSV"; then
  log_substep "Aggregating classification results"
  
  echo "Tool,Run,Fraction,SplitMode,eciRNA,EIciRNA,ciRNA,intergenic" > "$ALL_CLASS_CSV"
  
  find "$CLASS_DATA_DIR" -name "classification_results.csv" -type f | while read -r f; do
    tail -n +2 "$f" >> "$ALL_CLASS_CSV"
  done
  
  echo "  ✓ Aggregated $(( $(wc -l < "$ALL_CLASS_CSV") - 1 )) classification records"
fi

# Generate classification plots
if needs_run "${CLASS_PLOT_DIR}/gt_radial_f0.75_split.png"; then
  log_substep "Generating classification plots"
  
  cd "$CLASS_PLOT_DIR"
  cp "$ALL_CLASS_CSV" ./
  
  python3 - <<'PYPLOT'
import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

data = pd.read_csv('all_classification_results.csv')

CIRCRNA_COLORS = {
    'eciRNA': '#1b9e77',
    'EIciRNA': '#d95f02',
    'ciRNA': '#7570b3',
    'intergenic': '#e7298a'
}

TOOL_DISPLAY_NAMES = {
    'ciri-long': 'CIRI-long',
    'isocirc': 'IsoCirc',
    'circnick': 'CircNick-lrs',
    'ground_truth': 'Ground Truth'
}

def create_radial_bars(counts_dict, title, out_png, gt_counts_dict=None, gt_se_dict=None, tool_se_dict=None):
    TYPES = ["eciRNA","EIciRNA","ciRNA","intergenic"]
    vals = [int(counts_dict.get(t,0)) for t in TYPES]
    cols = [CIRCRNA_COLORS[t] for t in TYPES]
    
    is_gt = (gt_counts_dict is None)
    
    if is_gt:
        r_scaled = np.ones(len(TYPES))
        se_dict = gt_se_dict if gt_se_dict else {}
    else:
        gt_vals = [int(gt_counts_dict.get(t,0)) for t in TYPES]
        r_raw = np.array([v / gt if gt > 0 else 0 for v, gt in zip(vals, gt_vals)])
        r_scaled = np.sqrt(r_raw)
        se_dict = tool_se_dict if tool_se_dict else {}
    
    theta = np.linspace(0, 2*np.pi, len(TYPES), endpoint=False)
    width = 2*np.pi/len(TYPES) * 0.8
    
    fig = plt.figure(figsize=(14,14))
    ax = plt.subplot(111, polar=True)
    ax.set_theta_offset(np.pi/2)
    ax.set_theta_direction(-1)
    
    bars = ax.bar(theta, r_scaled, width=width, bottom=0, color=cols,
           edgecolor="black", linewidth=1.5, alpha=0.9)
    
    ax.set_xticks(theta)
    ax.set_xticklabels(TYPES, fontsize=18, fontweight="bold")
    
    if is_gt:
        ax.set_ylim(0, 1.15)
        ax.set_yticks([0, 0.5, 1.0])
        ax.set_yticklabels(['0%', '50%', '100%'], fontsize=14)
        
        for ang, typ, raw in zip(theta, TYPES, vals):
            if raw > 0:
                se = se_dict.get(typ, 0)
                
                if se > 0:
                    ax.plot([ang, ang], [0.48, 0.52], color='red', linewidth=3, zorder=10)
                    cap_width = 0.05
                    ax.plot([ang - cap_width, ang + cap_width], [0.48, 0.48],
                           color='red', linewidth=2, zorder=10)
                    ax.plot([ang - cap_width, ang + cap_width], [0.52, 0.52],
                           color='red', linewidth=2, zorder=10)
                
                label_text = f"{raw}"
                if se > 0:
                    label_text += f"\n±{se:.1f}"
                
                ax.text(ang, 0.5, label_text, ha="center", va="center",
                        fontsize=18, fontweight="bold", color="white",
                        bbox=dict(boxstyle='round,pad=0.5', facecolor='black', alpha=0.85),
                        zorder=11)
        
        ax.set_title(title + "\n(Each bar = 100% of that type)", 
                    fontsize=20, fontweight="bold", pad=25)
    else:
        ax.set_ylim(0, 1.15)
        tick_positions = np.sqrt([0, 0.1, 0.2, 0.4, 0.6, 0.8, 1.0])
        ax.set_yticks(tick_positions)
        ax.set_yticklabels(['0%', '10%', '20%', '40%', '60%', '80%', '100%'], fontsize=14)
        
        gt_vals = [int(gt_counts_dict.get(t,0)) for t in TYPES]
        r_raw = np.array([v / gt if gt > 0 else 0 for v, gt in zip(vals, gt_vals)])
        
        for ang, typ, raw, pct in zip(theta, TYPES, vals, r_raw):
            if pct > 0:
                se = se_dict.get(typ, 0)
                
                if se > 0:
                    ax.plot([ang, ang], [0.48, 0.52], color='red', linewidth=3, zorder=10)
                    cap_width = 0.05
                    ax.plot([ang - cap_width, ang + cap_width], [0.48, 0.48],
                           color='red', linewidth=2, zorder=10)
                    ax.plot([ang - cap_width, ang + cap_width], [0.52, 0.52],
                           color='red', linewidth=2, zorder=10)
                
                label_text = f"{raw}\n({pct*100:.0f}%)"
                if se > 0:
                    label_text = f"{raw}±{se:.1f}\n({pct*100:.0f}%)"
                
                ax.text(ang, 0.5, label_text, ha="center", va="center",
                       fontsize=16, fontweight="bold", color="white",
                       bbox=dict(boxstyle='round,pad=0.5', facecolor='black', alpha=0.85),
                       zorder=11)
        
        ax.set_title(title + "\n(% of GT detected per type)", 
                    fontsize=20, fontweight="bold", pad=25)
    
    if is_gt:
        fig.text(0.5, 0.02, 'Note: Each type scaled to 100%. Error bars show SE across runs.',
                ha="center", fontsize=13, style="italic", color='gray')
    else:
        fig.text(0.5, 0.02, 'Note: Sqrt scale. Each type shows % of GT detected. Error bars show SE.',
                ha="center", fontsize=13, style="italic", color='gray')
    
    plt.tight_layout(rect=[0, 0.04, 1, 1])
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close()

# Generate plots for each fraction/split combination
fractions = sorted(data['Fraction'].unique())
split_modes = sorted(data['SplitMode'].unique())

for fraction in fractions:
    for split_mode in split_modes:
        subset = data[(data['Fraction'] == fraction) & (data['SplitMode'] == split_mode)]
        if len(subset) == 0: continue
        
        # Ground truth
        gt_data = subset[subset['Tool'] == 'ground_truth']
        if len(gt_data) == 0: continue
        
        gt_avg = {t: gt_data[t].mean() for t in ['eciRNA', 'EIciRNA', 'ciRNA', 'intergenic']}
        gt_se = {t: gt_data[t].sem() if len(gt_data) > 1 else 0 for t in ['eciRNA', 'EIciRNA', 'ciRNA', 'intergenic']}
        
        title = f"Ground Truth - f={fraction} ({split_mode})"
        filename = f"gt_radial_f{fraction}_{split_mode}.png"
        create_radial_bars(gt_avg, title, filename, gt_se_dict=gt_se)
        
        # Tools
        for tool in ['ciri-long', 'isocirc', 'circnick']:
            tool_data = subset[subset['Tool'] == tool]
            if len(tool_data) == 0: continue
            
            tool_avg = {t: tool_data[t].mean() for t in ['eciRNA', 'EIciRNA', 'ciRNA', 'intergenic']}
            tool_se = {t: tool_data[t].sem() if len(tool_data) > 1 else 0 for t in ['eciRNA', 'EIciRNA', 'ciRNA', 'intergenic']}
            
            tool_display = TOOL_DISPLAY_NAMES.get(tool, tool)
            title = f"{tool_display} - f={fraction} ({split_mode})"
            filename = f"{tool}_radial_f{fraction}_{split_mode}.png"
            create_radial_bars(tool_avg, title, filename, gt_counts_dict=gt_avg, tool_se_dict=tool_se)

print("✓ Classification plots generated")
PYPLOT
  
  rm all_classification_results.csv
  cd - > /dev/null
  echo "  ✓ Classification plots complete"
fi

################################################################################
# STEP 4: TOOL PERFORMANCE BENCHMARKING
################################################################################
log_step "STEP 4: Tool Performance Benchmarking"

PERF_DATA_DIR="${OUTPUT_DIR}/data/performance"
PERF_PLOT_DIR="${OUTPUT_DIR}/plots/performance"
mkdir -p "$PERF_DATA_DIR" "$PERF_PLOT_DIR"

PERF_MARKER="${PERF_DATA_DIR}/.benchmarking_complete"

if [ ! -f "$PERF_MARKER" ]; then
  log_substep "Running tool benchmarking"
  
  benchmark_tool() {
    local run_name="$1"
    local tool_name="$2"
    local tool_bed="$3"
    local gt_bed="$4"
    local fraction="$5"
    local use_split="$6"
    
    local mode_dir="nosplit"
    [ "$use_split" = "true" ] && mode_dir="split"
    
    local out_dir="${PERF_DATA_DIR}/${run_name}/${mode_dir}/${tool_name}_f${fraction}"
    local metrics_file="${out_dir}/metrics.txt"
    
    if [ -f "$metrics_file" ]; then
      return 0
    fi
    
    mkdir -p "$out_dir"
    
    local SPLIT_OPT=""
    [ "$use_split" = "true" ] && SPLIT_OPT="-split"
    
    # Simple coord extraction
    grep -v '^#' "$tool_bed" 2>/dev/null | grep -v '^track' | grep -v '^browser' | \
      cut -f1-3 | sort -k1,1 -k2,2n -k3,3n | uniq > "${out_dir}/tool_coords.bed"
    
    grep -v '^#' "$gt_bed" 2>/dev/null | grep -v '^track' | grep -v '^browser' | \
      cut -f1-3 | sort -k1,1 -k2,2n -k3,3n | uniq > "${out_dir}/gt_coords.bed"
    
    local TOTAL_PRED=$(wc -l < "${out_dir}/tool_coords.bed" || echo 0)
    local GT_COUNT=$(wc -l < "${out_dir}/gt_coords.bed" || echo 0)
    
    if [ "$TOTAL_PRED" -eq 0 ]; then
      TP=0; FP=0; FN="$GT_COUNT"
    else
      bedtools intersect -a "${out_dir}/tool_coords.bed" -b "${out_dir}/gt_coords.bed" \
        -f "$fraction" -r -wa -u ${SPLIT_OPT} > "${out_dir}/TP.bed" 2>/dev/null || touch "${out_dir}/TP.bed"
      bedtools intersect -a "${out_dir}/tool_coords.bed" -b "${out_dir}/gt_coords.bed" \
        -f "$fraction" -r -v ${SPLIT_OPT} > "${out_dir}/FP.bed" 2>/dev/null || touch "${out_dir}/FP.bed"
      bedtools intersect -a "${out_dir}/gt_coords.bed" -b "${out_dir}/tool_coords.bed" \
        -f "$fraction" -r -v ${SPLIT_OPT} > "${out_dir}/FN.bed" 2>/dev/null || touch "${out_dir}/FN.bed"
      
      local TP=$(wc -l < "${out_dir}/TP.bed")
      local FP=$(wc -l < "${out_dir}/FP.bed")
      local FN=$(wc -l < "${out_dir}/FN.bed")
    fi
    
    local PREC REC F1
    if [ $((TP + FP)) -ne 0 ]; then
      PREC=$(echo "scale=6; ${TP} / (${TP} + ${FP})" | bc)
    else
      PREC="0"
    fi
    if [ $((TP + FN)) -ne 0 ]; then
      REC=$(echo "scale=6; ${TP} / (${TP} + ${FN})" | bc)
    else
      REC="0"
    fi
    if [ "$(echo "${PREC} + ${REC} > 0" | bc)" -eq 1 ]; then
      F1=$(echo "scale=6; 2 * (${PREC} * ${REC}) / (${PREC} + ${REC})" | bc)
    else
      F1="0"
    fi
    
    cat > "$metrics_file" <<EOF
Run: ${run_name}
Tool: ${tool_name}
Fraction: ${fraction}
Split: ${use_split}
TP: ${TP}
FP: ${FP}
FN: ${FN}
Precision: ${PREC}
Recall: ${REC}
F1_Score: ${F1}
EOF
  }
  
  for rd in $(run_dirs); do
    run_name="$(basename "$rd")"
    sanitized_run="${SANITIZED_DIR}/${run_name}"
    
    gt="${sanitized_run}/${GT_BED_NAME}"
    [ -f "$gt" ] || continue
    
    for tool_name in ciri-long isocirc circnick; do
      case $tool_name in
        ciri-long) tool_bed="${sanitized_run}/${CIRI_BED_NAME}" ;;
        isocirc) tool_bed="${sanitized_run}/${ISO_BED_NAME}" ;;
        circnick) tool_bed="${sanitized_run}/${NICK_BED_NAME}" ;;
      esac
      
      [ -f "$tool_bed" ] || continue
      
      for frac in "${FRACTIONS[@]}"; do
        for split in "${SPLIT_MODES[@]}"; do
          benchmark_tool "$run_name" "$tool_name" "$tool_bed" "$gt" "$frac" "$split"
        done
      done
    done
  done
  
  touch "$PERF_MARKER"
  echo "  ✓ Tool benchmarking complete"
fi

# Aggregate performance metrics
ALL_PERF_CSV="${PERF_DATA_DIR}/all_tool_metrics.csv"

if needs_run "$ALL_PERF_CSV"; then
  log_substep "Aggregating performance metrics"
  
  echo "Run,Tool,Fraction,Split,TP,FP,FN,Precision,Recall,F1_Score" > "$ALL_PERF_CSV"
  
  find "$PERF_DATA_DIR" -name "metrics.txt" -type f | while read -r mfile; do
    RUN=$(grep "^Run:" "$mfile" | awk '{print $2}')
    TOOL=$(grep "^Tool:" "$mfile" | awk '{print $2}')
    FRAC=$(grep "^Fraction:" "$mfile" | awk '{print $2}')
    SPLIT=$(grep "^Split:" "$mfile" | awk '{print $2}')
    TP=$(grep "^TP:" "$mfile" | awk '{print $2}')
    FP=$(grep "^FP:" "$mfile" | awk '{print $2}')
    FN=$(grep "^FN:" "$mfile" | awk '{print $2}')
    PREC=$(grep "^Precision:" "$mfile" | awk '{print $2}')
    REC=$(grep "^Recall:" "$mfile" | awk '{print $2}')
    F1=$(grep "^F1_Score:" "$mfile" | awk '{print $2}')
    
    [ -n "$RUN" ] && [ -n "$TOOL" ] && \
      echo "$RUN,$TOOL,$FRAC,$SPLIT,$TP,$FP,$FN,$PREC,$REC,$F1" >> "$ALL_PERF_CSV"
  done
  
  echo "  ✓ Aggregated $(( $(wc -l < "$ALL_PERF_CSV") - 1 )) performance records"
fi

# Generate performance plots (main + heatmap with new color scheme)
if needs_run "${PERF_PLOT_DIR}/tool_performance.png"; then
  log_substep "Generating performance plots"
  
  cd "$PERF_PLOT_DIR"
  cp "$ALL_PERF_CSV" ./
  
  Rscript - <<'RSCRIPT'
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(scales)
})

data <- read.csv("all_tool_metrics.csv", stringsAsFactors=FALSE)

data <- data %>%
  mutate(
    Fraction = factor(Fraction, levels=c("0.25","0.5","0.75","0.95")),
    Split = as.logical(Split),
    ToolName = case_when(
      Tool == "ciri-long" ~ "CIRI-long",
      Tool == "isocirc"   ~ "IsoCirc",
      Tool == "circnick"  ~ "CircNick-lrs",
      TRUE ~ Tool
    ),
    SplitMode = ifelse(Split, "exon-level", "transcriptome-level")
  )

tool_colors <- c("CIRI-long"="#0072B2","IsoCirc"="#009E73","CircNick-lrs"="#D55E00")
data$ToolName <- factor(data$ToolName, levels=names(tool_colors))

summary_data <- data %>%
  group_by(Tool, ToolName, Fraction, Split, SplitMode) %>%
  summarize(
    Precision_Mean = mean(Precision, na.rm=TRUE),
    Recall_Mean    = mean(Recall,    na.rm=TRUE),
    F1_Score_Mean  = mean(F1_Score,  na.rm=TRUE),
    Precision_SE   = sd(Precision, na.rm=TRUE)/sqrt(n()),
    Recall_SE      = sd(Recall,    na.rm=TRUE)/sqrt(n()),
    F1_Score_SE    = sd(F1_Score,  na.rm=TRUE)/sqrt(n()),
    .groups="drop"
  )

plot_data <- summary_data %>%
  pivot_longer(cols=c(Precision_Mean, Recall_Mean, F1_Score_Mean),
               names_to="Metric", values_to="Value") %>%
  mutate(Metric_Base=gsub("_Mean","",Metric)) %>%
  left_join(
    summary_data %>%
      pivot_longer(cols=c(Precision_SE, Recall_SE, F1_Score_SE),
                   names_to="SE_Metric", values_to="SE") %>%
      mutate(Metric_Base=gsub("_SE","",SE_Metric)),
    by=c("ToolName","Fraction","Split","SplitMode","Metric_Base")
  )

plot_data$Metric_Base <- factor(plot_data$Metric_Base, levels=c("Precision","Recall","F1_Score"))
plot_data$ShapeKey <- interaction(plot_data$SplitMode, plot_data$Fraction, sep=".")

shape_vals <- c(
  "exon-level.0.25"=15, "exon-level.0.5"=16, "exon-level.0.75"=17, "exon-level.0.95"=18,
  "transcriptome-level.0.25"=0, "transcriptome-level.0.5"=1, 
  "transcriptome-level.0.75"=2, "transcriptome-level.0.95"=5
)

png("tool_performance.png", width=1800, height=1200, res=150)
p <- ggplot(plot_data, aes(x=Metric_Base, y=Value, color=ToolName, shape=ShapeKey)) +
  geom_errorbar(aes(ymin=Value-SE, ymax=Value+SE), width=0.2, position=position_dodge(width=0.6)) +
  geom_point(size=4, stroke=1.2, position=position_dodge(width=0.6)) +
  scale_color_manual(values=tool_colors) +
  scale_shape_manual(values=shape_vals, name="Overlap fraction\n(filled=exon-level)") +
  scale_y_continuous(limits=c(0,1), labels=scales::percent_format(accuracy=1)) +
  labs(title="Tool Performance Metrics (Mean ± SE)", x="Metric", y="Value", color="Tool") +
  theme_minimal(base_size=14) +
  theme(plot.title=element_text(size=18, face="bold", hjust=0.5),
        legend.position="bottom")
print(p); dev.off()

# Contingency heatmap with yellow-to-blue color scheme
cont_sum <- data %>%
  group_by(ToolName, Fraction, SplitMode) %>%
  summarize(TP=sum(TP,na.rm=TRUE), FP=sum(FP,na.rm=TRUE), FN=sum(FN,na.rm=TRUE), .groups="drop")

heat <- cont_sum %>% pivot_longer(cols=c(TP,FP,FN), names_to="Metric", values_to="Count")
heat$Metric <- factor(heat$Metric, levels=c("TP","FP","FN"))
heat$ToolName <- factor(heat$ToolName, levels=names(tool_colors))

png("tool_contingency_heatmap.png", width=1800, height=1200, res=150)
pH <- ggplot(heat, aes(x=Fraction, y=ToolName, fill=Count)) +
  geom_tile(color="white") +
  facet_grid(SplitMode ~ Metric) +
  scale_fill_gradient(low="#FFFFCC", high="#08519C", trans="log10",
                      breaks=c(10,100,1000,10000), labels=scales::comma) +
  geom_text(aes(label=scales::comma(Count)), size=3) +
  labs(title="Tool Detection Performance: Contingency Table",
       subtitle="TP/FP/FN totals across runs",
       x="Overlap Fraction", y="Tool", fill="Count (log)") +
  theme_minimal() +
  theme(plot.title=element_text(size=14, face="bold", hjust=0.5),
        plot.subtitle=element_text(size=11, hjust=0.5, face="italic"),
        strip.text=element_text(size=10, face="bold"))
print(pH); dev.off()

cat("✓ Performance plots generated\n")
RSCRIPT
  
  rm all_tool_metrics.csv
  cd - > /dev/null
  echo "  ✓ Performance plots complete"
fi

################################################################################
# STEP 5: COMBO PERFORMANCE
################################################################################
log_step "STEP 5: Combo Performance"

COMBO_DATA_DIR="${OUTPUT_DIR}/data/combos"
COMBO_PLOT_DIR="${OUTPUT_DIR}/plots/combos"
mkdir -p "$COMBO_DATA_DIR" "$COMBO_PLOT_DIR"

COMBO_MARKER="${COMBO_DATA_DIR}/.combos_complete"

if [ ! -f "$COMBO_MARKER" ]; then
  log_substep "Preparing and benchmarking combos"
  
  for rd in $(run_dirs); do
    run_name="$(basename "$rd")"
    sanitized_run="${SANITIZED_DIR}/${run_name}"
    
    gt="${sanitized_run}/${GT_BED_NAME}"
    [ -f "$gt" ] || continue
    
    # Prepare combo beds
    run_combo_dir="${COMBO_DATA_DIR}/${run_name}"
    mkdir -p "${run_combo_dir}/combos"
    
    ciri="${sanitized_run}/${CIRI_BED_NAME}"
    iso="${sanitized_run}/${ISO_BED_NAME}"
    nick="${sanitized_run}/${NICK_BED_NAME}"
    
    # Copy GT
    cp "$gt" "${run_combo_dir}/gt.bed"
    
    # Create combos (intersections)
    if [ -f "$ciri" ] && [ -f "$iso" ]; then
      bedtools intersect -a "$ciri" -b "$iso" -wa 2>/dev/null | sort | uniq \
        > "${run_combo_dir}/combos/ciri_iso.bed"
    else
      touch "${run_combo_dir}/combos/ciri_iso.bed"
    fi
    
    if [ -f "$ciri" ] && [ -f "$nick" ]; then
      bedtools intersect -a "$ciri" -b "$nick" -wa 2>/dev/null | sort | uniq \
        > "${run_combo_dir}/combos/ciri_nick.bed"
    else
      touch "${run_combo_dir}/combos/ciri_nick.bed"
    fi
    
    if [ -f "$iso" ] && [ -f "$nick" ]; then
      bedtools intersect -a "$iso" -b "$nick" -wa 2>/dev/null | sort | uniq \
        > "${run_combo_dir}/combos/iso_nick.bed"
    else
      touch "${run_combo_dir}/combos/iso_nick.bed"
    fi
    
    # Union
    cat /dev/null > "${run_combo_dir}/combos/all_union.bed"
    [ -f "$ciri" ] && cat "$ciri" >> "${run_combo_dir}/combos/all_union.bed"
    [ -f "$iso" ] && cat "$iso" >> "${run_combo_dir}/combos/all_union.bed"
    [ -f "$nick" ] && cat "$nick" >> "${run_combo_dir}/combos/all_union.bed"
    sort "${run_combo_dir}/combos/all_union.bed" | uniq > "${run_combo_dir}/combos/all_union.bed.tmp"
    mv "${run_combo_dir}/combos/all_union.bed.tmp" "${run_combo_dir}/combos/all_union.bed"
    
    # Benchmark each combo
    for combo in ciri_iso ciri_nick iso_nick all_union; do
      combo_bed="${run_combo_dir}/combos/${combo}.bed"
      [ -s "$combo_bed" ] || continue
      
      for frac in "${FRACTIONS[@]}"; do
        for split in "${SPLIT_MODES[@]}"; do
          mode_dir="nosplit"
          [ "$split" = "true" ] && mode_dir="split"
          
          out_dir="${run_combo_dir}/${mode_dir}/${combo}_f${frac}"
          metrics_file="${out_dir}/metrics.txt"
          
          [ -f "$metrics_file" ] && continue
          
          mkdir -p "$out_dir"
          
          SPLIT_OPT=""
          [ "$split" = "true" ] && SPLIT_OPT="-split"
          
          # Extract coords
          grep -v '^#' "$combo_bed" | grep -v '^track' | grep -v '^browser' | \
            cut -f1-3 | sort -k1,1 -k2,2n -k3,3n | uniq > "${out_dir}/combo_coords.bed"
          
          grep -v '^#' "$gt" | grep -v '^track' | grep -v '^browser' | \
            cut -f1-3 | sort -k1,1 -k2,2n -k3,3n | uniq > "${out_dir}/gt_coords.bed"
          
          TOTAL_PRED=$(wc -l < "${out_dir}/combo_coords.bed")
          GT_COUNT=$(wc -l < "${out_dir}/gt_coords.bed")
          
          if [ "$TOTAL_PRED" -eq 0 ]; then
            TP=0; FP=0; FN="$GT_COUNT"
          else
            bedtools intersect -a "${out_dir}/combo_coords.bed" -b "${out_dir}/gt_coords.bed" \
              -f "$frac" -r -wa -u ${SPLIT_OPT} > "${out_dir}/TP.bed" 2>/dev/null || touch "${out_dir}/TP.bed"
            bedtools intersect -a "${out_dir}/combo_coords.bed" -b "${out_dir}/gt_coords.bed" \
              -f "$frac" -r -v ${SPLIT_OPT} > "${out_dir}/FP.bed" 2>/dev/null || touch "${out_dir}/FP.bed"
            bedtools intersect -a "${out_dir}/gt_coords.bed" -b "${out_dir}/combo_coords.bed" \
              -f "$frac" -r -v ${SPLIT_OPT} > "${out_dir}/FN.bed" 2>/dev/null || touch "${out_dir}/FN.bed"
            
            TP=$(wc -l < "${out_dir}/TP.bed")
            FP=$(wc -l < "${out_dir}/FP.bed")
            FN=$(wc -l < "${out_dir}/FN.bed")
          fi
          
          PREC="0"; REC="0"; F1="0"
          if [ $((TP + FP)) -ne 0 ]; then
            PREC=$(echo "scale=6; ${TP} / (${TP} + ${FP})" | bc)
          fi
          if [ $((TP + FN)) -ne 0 ]; then
            REC=$(echo "scale=6; ${TP} / (${TP} + ${FN})" | bc)
          fi
          if [ "$(echo "${PREC} + ${REC} > 0" | bc)" -eq 1 ]; then
            F1=$(echo "scale=6; 2 * (${PREC} * ${REC}) / (${PREC} + ${REC})" | bc)
          fi
          
          cat > "$metrics_file" <<EOF
Run: ${run_name}
Combo: ${combo}
Fraction: ${frac}
Split: ${split}
TP: ${TP}
FP: ${FP}
FN: ${FN}
Precision: ${PREC}
Recall: ${REC}
F1_Score: ${F1}
EOF
        done
      done
    done
  done
  
  touch "$COMBO_MARKER"
  echo "  ✓ Combo benchmarking complete"
fi

# Aggregate combo metrics
ALL_COMBO_CSV="${COMBO_DATA_DIR}/all_combo_metrics.csv"

if needs_run "$ALL_COMBO_CSV"; then
  log_substep "Aggregating combo metrics"
  
  echo "Run,Combo,Fraction,Split,TP,FP,FN,Precision,Recall,F1_Score" > "$ALL_COMBO_CSV"
  
  find "$COMBO_DATA_DIR" -name "metrics.txt" -type f | while read -r mfile; do
    RUN=$(grep "^Run:" "$mfile" | awk '{print $2}')
    COMBO=$(grep "^Combo:" "$mfile" | awk '{print $2}')
    FRAC=$(grep "^Fraction:" "$mfile" | awk '{print $2}')
    SPLIT=$(grep "^Split:" "$mfile" | awk '{print $2}')
    TP=$(grep "^TP:" "$mfile" | awk '{print $2}')
    FP=$(grep "^FP:" "$mfile" | awk '{print $2}')
    FN=$(grep "^FN:" "$mfile" | awk '{print $2}')
    PREC=$(grep "^Precision:" "$mfile" | awk '{print $2}')
    REC=$(grep "^Recall:" "$mfile" | awk '{print $2}')
    F1=$(grep "^F1_Score:" "$mfile" | awk '{print $2}')
    
    [ -n "$RUN" ] && [ -n "$COMBO" ] && \
      echo "$RUN,$COMBO,$FRAC,$SPLIT,$TP,$FP,$FN,$PREC,$REC,$F1" >> "$ALL_COMBO_CSV"
  done
  
  echo "  ✓ Aggregated $(( $(wc -l < "$ALL_COMBO_CSV") - 1 )) combo records"
fi

# Generate combo plot (combinations_performance.png only, with fixed legend)
if needs_run "${COMBO_PLOT_DIR}/combinations_performance.png"; then
  log_substep "Generating combo plot"
  
  cd "$COMBO_PLOT_DIR"
  cp "$ALL_COMBO_CSV" ./
  
  Rscript - <<'RSCRIPT'
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(scales)
})

data <- read.csv("all_combo_metrics.csv", stringsAsFactors=FALSE)

data <- data %>%
  mutate(
    Fraction = factor(Fraction, levels=c("0.25","0.5","0.75","0.95")),
    Split = as.logical(Split),
    ComboName = case_when(
      Combo == "ciri_iso"  ~ "CIRI-long ∩ IsoCirc",
      Combo == "ciri_nick" ~ "CIRI-long ∩ CircNick",
      Combo == "iso_nick"  ~ "IsoCirc ∩ CircNick",
      Combo == "all_union" ~ "Union (All Tools)",
      TRUE ~ Combo
    ),
    SplitMode = ifelse(Split, "exon-level", "transcriptome-level")
  )

combo_colors <- c(
  "CIRI-long ∩ IsoCirc"   = "#4575B4",
  "CIRI-long ∩ CircNick"  = "#91BFDB",
  "IsoCirc ∩ CircNick"    = "#FC8D59",
  "Union (All Tools)"     = "#D73027"
)

data$ComboName <- factor(data$ComboName, levels=names(combo_colors))

summary_data <- data %>%
  group_by(Combo, ComboName, Fraction, Split, SplitMode) %>%
  summarize(
    Precision_Mean = mean(Precision, na.rm=TRUE),
    Recall_Mean    = mean(Recall,    na.rm=TRUE),
    F1_Score_Mean  = mean(F1_Score,  na.rm=TRUE),
    Precision_SE   = sd(Precision, na.rm=TRUE)/sqrt(n()),
    Recall_SE      = sd(Recall,    na.rm=TRUE)/sqrt(n()),
    F1_Score_SE    = sd(F1_Score,  na.rm=TRUE)/sqrt(n()),
    .groups="drop"
  )

plot_data <- summary_data %>%
  pivot_longer(cols=c(Precision_Mean, Recall_Mean, F1_Score_Mean),
               names_to="Metric", values_to="Value") %>%
  mutate(Metric_Base=gsub("_Mean","",Metric)) %>%
  left_join(
    summary_data %>%
      pivot_longer(cols=c(Precision_SE, Recall_SE, F1_Score_SE),
                   names_to="SE_Metric", values_to="SE") %>%
      mutate(Metric_Base=gsub("_SE","",SE_Metric)),
    by=c("ComboName","Fraction","Split","SplitMode","Metric_Base")
  )

plot_data$Metric_Base <- factor(plot_data$Metric_Base, levels=c("Precision","Recall","F1_Score"))
plot_data$ShapeKey <- interaction(plot_data$SplitMode, plot_data$Fraction, sep=".")

shape_vals <- c(
  "exon-level.0.25"=15, "exon-level.0.5"=16, "exon-level.0.75"=17, "exon-level.0.95"=18,
  "transcriptome-level.0.25"=0, "transcriptome-level.0.5"=1, 
  "transcriptome-level.0.75"=2, "transcriptome-level.0.95"=5
)

png("combinations_performance.png", width=2000, height=1200, res=150)
p <- ggplot(plot_data, aes(x=Metric_Base, y=Value, color=ComboName, shape=ShapeKey)) +
  geom_errorbar(aes(ymin=Value-SE, ymax=Value+SE), width=0.2, position=position_dodge(width=0.6)) +
  geom_point(size=4, stroke=1.2, position=position_dodge(width=0.6)) +
  scale_color_manual(values=combo_colors, name="Combination") +
  scale_shape_manual(values=shape_vals, name="Overlap fraction\n(filled=exon-level)") +
  scale_y_continuous(limits=c(0,1), labels=scales::percent_format(accuracy=1)) +
  labs(title="Combo Performance Metrics (Mean ± SE)", x="Metric", y="Value") +
  theme_minimal(base_size=14) +
  theme(plot.title=element_text(size=18, face="bold", hjust=0.5),
        legend.position="bottom",
        legend.box="vertical",
        legend.margin=margin(t=10))
print(p); dev.off()

cat("✓ Combo plot generated\n")
RSCRIPT
  
  rm all_combo_metrics.csv
  cd - > /dev/null
  echo "  ✓ Combo plot complete"
fi

################################################################################
# STEP 6: EXPRESSION ANALYSIS
################################################################################
log_step "STEP 6: Expression Analysis"

EXPR_DATA_DIR="${OUTPUT_DIR}/data/expression"
EXPR_PLOT_DIR="${OUTPUT_DIR}/plots/expression"
mkdir -p "$EXPR_DATA_DIR" "$EXPR_PLOT_DIR"

EXPR_MARKER="${EXPR_DATA_DIR}/.expression_complete"

if [ ! -f "$EXPR_MARKER" ]; then
  log_substep "Parsing expression data"
  
  # Parse expression for each run
  for rd in $(run_dirs); do
    run_name="$(basename "$rd")"
    sanitized_run="${SANITIZED_DIR}/${run_name}"
    
    out_expr="${EXPR_DATA_DIR}/${run_name}"
    mkdir -p "$out_expr"
    
    # Ground truth - convert BED (0-based) to 1-based
    gt_bed="${sanitized_run}/${GT_BED_NAME}"
    if [ -f "$gt_bed" ]; then
      echo "coord,expression,tool,run" > "${out_expr}/gt_expression.csv"
      awk 'BEGIN{OFS=","} !/^#/ && !/^track/ && !/^browser/ {
        coord=$1":"($2+1)"-"$3; expr=$5; if(expr=="") expr=0; 
        print coord, expr, "ground_truth", "'"$run_name"'"
      }' "$gt_bed" >> "${out_expr}/gt_expression.csv"
    fi
    
    # CIRI-long (already 1-based)
    ciri_expr="${sanitized_run}/${CIRI_EXPR_NAME}"
    if [ -f "$ciri_expr" ]; then
      echo "coord,expression,tool,run" > "${out_expr}/ciri_expression.csv"
      awk 'BEGIN{FS="\t"; OFS=","} NR>1 {
        coord=$1; expr=$2; if(expr=="") expr=0; 
        print coord, expr, "ciri-long", "'"$run_name"'"
      }' "$ciri_expr" >> "${out_expr}/ciri_expression.csv"
    fi
    
    # IsoCirc - convert to 1-based
    iso_expr="${sanitized_run}/${ISO_EXPR_NAME}"
    if [ -f "$iso_expr" ]; then
      echo "coord,expression,tool,run" > "${out_expr}/iso_expression.csv"
      awk 'BEGIN{FS="\t"; OFS=","} !/^#/ {
        coord=$2":"($3+1)"-"$4; expr=$(NF-1); if(expr=="") expr=0; 
        print coord, expr, "isocirc", "'"$run_name"'"
      }' "$iso_expr" >> "${out_expr}/iso_expression.csv"
    fi
    
    # CircNick - convert to 1-based
    nick_expr="${sanitized_run}/${NICK_EXPR_NAME}"
    if [ -f "$nick_expr" ]; then
      echo "coord,expression,tool,run" > "${out_expr}/circnick_expression.csv"
      awk 'BEGIN{FS="\t"; OFS=","} NR>1 {
        coord=$2":"($3+1)"-"$4; expr=$6; if(expr=="") expr=0; 
        print coord, expr, "circnick", "'"$run_name"'"
      }' "$nick_expr" >> "${out_expr}/circnick_expression.csv"
    fi
  done
  
  touch "$EXPR_MARKER"
  echo "  ✓ Expression data parsed"
fi

# Aggregate expression data - IMPROVED with direct concatenation
ALL_EXPR_CSV="${EXPR_DATA_DIR}/all_expression.csv"

if needs_run "$ALL_EXPR_CSV"; then
  log_substep "Aggregating expression data"
  
  echo "coord,expression,tool,run" > "$ALL_EXPR_CSV"
  
  # Direct concatenation - much faster than find
  for rd in $(run_dirs); do
    run_name="$(basename "$rd")"
    run_expr="${EXPR_DATA_DIR}/${run_name}"
    
    if [ -d "$run_expr" ]; then
      echo "    Processing ${run_name}..."
      
      # Concatenate each tool's expression file
      for expr_file in "${run_expr}/gt_expression.csv" \
                       "${run_expr}/ciri_expression.csv" \
                       "${run_expr}/iso_expression.csv" \
                       "${run_expr}/circnick_expression.csv"; do
        if [ -f "$expr_file" ]; then
          tail -n +2 "$expr_file" >> "$ALL_EXPR_CSV"
        fi
      done
    fi
  done
  
  record_count=$(( $(wc -l < "$ALL_EXPR_CSV") - 1 ))
  echo "  ✓ Aggregated ${record_count} expression records"
fi

# Generate expression plots
if needs_run "${EXPR_PLOT_DIR}/expression_correlation.png"; then
  log_substep "Generating expression plots"
  
  cd "$EXPR_PLOT_DIR"
  cp "$ALL_EXPR_CSV" ./
  
  python3 - <<'PYEXPR'
import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats

print("Loading expression data...")
data = pd.read_csv('all_expression.csv')
print(f"Loaded {len(data)} records")

# Merge GT with each tool
print("Merging with ground truth...")
gt_data = data[data['tool'] == 'ground_truth'][['coord', 'expression', 'run']].rename(columns={'expression': 'gt_expr'})

tool_names = {'ciri-long': 'CIRI-long', 'isocirc': 'IsoCirc', 'circnick': 'CircNick-lrs'}
tool_colors = {'ciri-long': '#0072B2', 'isocirc': '#009E73', 'circnick': '#D55E00'}

fig, axes = plt.subplots(1, 3, figsize=(18, 6))
fig.suptitle('Expression Correlation: Tool vs Ground Truth', fontsize=18, fontweight='bold')

for idx, (tool, tool_label) in enumerate(tool_names.items()):
    print(f"Processing {tool}...")
    ax = axes[idx]
    
    tool_data = data[data['tool'] == tool][['coord', 'expression', 'run']].rename(columns={'expression': f'{tool}_expr'})
    
    merged = pd.merge(gt_data, tool_data, on=['coord', 'run'], how='left')
    merged[f'{tool}_expr'] = merged[f'{tool}_expr'].fillna(0)
    
    # Filter non-zero in both
    nonzero = merged[(merged['gt_expr'] > 0) & (merged[f'{tool}_expr'] > 0)]
    
    if len(nonzero) > 10:
        ax.scatter(nonzero['gt_expr'], nonzero[f'{tool}_expr'], alpha=0.4, s=20, 
                  color=tool_colors[tool], edgecolors='black', linewidth=0.2)
        
        # Perfect line
        max_val = max(nonzero['gt_expr'].max(), nonzero[f'{tool}_expr'].max())
        ax.plot([0, max_val], [0, max_val], 'r--', lw=2, alpha=0.7, label='Perfect (y=x)')
        
        # Correlation
        r, p = stats.pearsonr(nonzero['gt_expr'], nonzero[f'{tool}_expr'])
        r2 = r ** 2
        
        ax.text(0.05, 0.95, f"Pearson r = {r:.3f}\nR² = {r2:.3f}\nn = {len(nonzero):,}", 
               transform=ax.transAxes, va='top', fontsize=11, fontweight='bold',
               bbox=dict(boxstyle='round', facecolor='white', alpha=0.9, edgecolor='black'))
        
        ax.legend(loc='lower right', fontsize=10)
        ax.set_xlabel('Ground Truth Expression', fontsize=12, fontweight='bold')
        ax.set_ylabel(f'{tool_label} Predicted', fontsize=12, fontweight='bold')
        ax.set_title(f'{tool_label}', fontsize=14, fontweight='bold')
        ax.grid(True, alpha=0.3, linestyle='--')
        ax.set_xlim(left=0); ax.set_ylim(bottom=0)

print("Saving plot...")
plt.tight_layout()
plt.savefig('expression_correlation.png', dpi=300, bbox_inches='tight')
plt.close()

print("✓ Expression plots generated")
PYEXPR
  
  rm all_expression.csv
  cd - > /dev/null
  echo "  ✓ Expression plots complete"
fi

################################################################################
# STEP 7: UPSET ANALYSIS
################################################################################
log_step "STEP 7: UpSet Analysis"

UPSET_DATA_DIR="${OUTPUT_DIR}/data/upset"
UPSET_PLOT_DIR="${OUTPUT_DIR}/plots/upset"
mkdir -p "$UPSET_DATA_DIR" "$UPSET_PLOT_DIR"

# Use fraction 0.75 for UpSet analysis
UPSET_FRAC="0.75"

UPSET_MARKER="${UPSET_DATA_DIR}/.upset_complete"

if [ ! -f "$UPSET_MARKER" ]; then
  log_substep "Preparing UpSet data (TP and all predictions)"
  
  for rd in $(run_dirs); do
    run_name="$(basename "$rd")"
    sanitized_run="${SANITIZED_DIR}/${run_name}"
    
    out_upset="${UPSET_DATA_DIR}/${run_name}"
    mkdir -p "$out_upset"
    
    gt="${sanitized_run}/${GT_BED_NAME}"
    [ -f "$gt" ] || continue
    
    # Get GT count
    grep -v '^#' "$gt" | grep -v '^track' | grep -v '^browser' | cut -f1-3 | sort -u | wc -l > "${out_upset}/gt_count.txt"
    
    # Get TP for each tool from performance data
    for tool in ciri-long isocirc circnick; do
      tp_file="${PERF_DATA_DIR}/${run_name}/nosplit/${tool}_f${UPSET_FRAC}/TP.bed"
      
      if [ -f "$tp_file" ]; then
        awk '{print $1":"($2+1)"-"$3}' "$tp_file" | sort -u > "${out_upset}/${tool}_tp_ids.txt"
      else
        touch "${out_upset}/${tool}_tp_ids.txt"
      fi
    done
    
    # Get all predictions for each tool
    for tool in ciri-long isocirc circnick; do
      case $tool in
        ciri-long) tool_bed="${sanitized_run}/${CIRI_BED_NAME}" ;;
        isocirc) tool_bed="${sanitized_run}/${ISO_BED_NAME}" ;;
        circnick) tool_bed="${sanitized_run}/${NICK_BED_NAME}" ;;
      esac
      
      if [ -f "$tool_bed" ]; then
        grep -v '^#' "$tool_bed" | grep -v '^track' | grep -v '^browser' | \
          awk '{print $1":"($2+1)"-"$3}' | sort -u > "${out_upset}/${tool}_all_ids.txt"
      else
        touch "${out_upset}/${tool}_all_ids.txt"
      fi
    done
  done
  
  touch "$UPSET_MARKER"
  echo "  ✓ UpSet data prepared"
fi

# Generate UpSet plots
if needs_run "${UPSET_PLOT_DIR}/combined_upset_plot.png"; then
  log_substep "Generating UpSet plot"
  
  cd "$UPSET_PLOT_DIR"
  
  # Combine data from all runs
  python3 - <<'PYUPSET'
import os
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import patches
import numpy as np

data_dir = os.environ.get('UPSET_DATA_DIR', '../data/upset')

# Collect TP IDs and all IDs from all runs
tp_ids = {'ciri-long': set(), 'isocirc': set(), 'circnick': set()}
all_ids = {'ciri-long': set(), 'isocirc': set(), 'circnick': set()}
gt_total = 0

for run_dir in os.listdir(data_dir):
    run_path = os.path.join(data_dir, run_dir)
    if not os.path.isdir(run_path): continue
    
    # GT count
    gt_file = os.path.join(run_path, 'gt_count.txt')
    if os.path.exists(gt_file):
        with open(gt_file) as f:
            gt_total = max(gt_total, int(f.read().strip()))
    
    # TP IDs
    for tool in ['ciri-long', 'isocirc', 'circnick']:
        tp_file = os.path.join(run_path, f'{tool}_tp_ids.txt')
        if os.path.exists(tp_file):
            with open(tp_file) as f:
                tp_ids[tool].update(line.strip() for line in f if line.strip())
        
        all_file = os.path.join(run_path, f'{tool}_all_ids.txt')
        if os.path.exists(all_file):
            with open(all_file) as f:
                all_ids[tool].update(line.strip() for line in f if line.strip())

print(f"GT total: {gt_total}")
print(f"TP counts: CIRI={len(tp_ids['ciri-long'])}, Iso={len(tp_ids['isocirc'])}, Nick={len(tp_ids['circnick'])}")
print(f"All counts: CIRI={len(all_ids['ciri-long'])}, Iso={len(all_ids['isocirc'])}, Nick={len(all_ids['circnick'])}")

# Create UpSetR-style plot for TP
def create_upset_plot(id_dict, title, output_file, gt_total):
    tools = ['CIRI-long', 'IsoCirc', 'CircNick-lrs']
    tool_keys = ['ciri-long', 'isocirc', 'circnick']
    
    # Get all unique IDs
    all_unique = set()
    for ids in id_dict.values():
        all_unique.update(ids)
    
    # Create presence/absence matrix
    matrix = {}
    for coord in all_unique:
        presence = tuple(1 if coord in id_dict[key] else 0 for key in tool_keys)
        if presence not in matrix:
            matrix[presence] = 0
        matrix[presence] += 1
    
    # Sort by count descending
    sorted_combos = sorted(matrix.items(), key=lambda x: x[1], reverse=True)
    
    # Prepare data for plotting
    combinations = [combo for combo, count in sorted_combos]
    counts = [count for combo, count in sorted_combos]
    percentages = [100 * count / gt_total for count in counts]
    
    # Set sizes
    n_tools = len(tools)
    n_combos = len(combinations)
    
    fig = plt.figure(figsize=(16, 10))
    
    # Create grid
    gs = fig.add_gridspec(3, 2, height_ratios=[3, 0.5, 1.5], width_ratios=[1, 6],
                          hspace=0.05, wspace=0.05)
    
    # Main bar chart (top right)
    ax_bars = fig.add_subplot(gs[0, 1])
    x_pos = np.arange(n_combos)
    bars = ax_bars.bar(x_pos, counts, color='#313695', edgecolor='black', linewidth=0.8)
    
    # Add count labels on bars
    for i, (count, pct) in enumerate(zip(counts, percentages)):
        ax_bars.text(i, count + max(counts)*0.02, f'{count}\n({pct:.1f}%)', 
                    ha='center', va='bottom', fontsize=9, fontweight='bold')
    
    ax_bars.set_ylabel('Intersection Size', fontsize=14, fontweight='bold')
    ax_bars.set_ylim(0, max(counts) * 1.15)
    ax_bars.set_xlim(-0.5, n_combos - 0.5)
    ax_bars.set_xticks([])
    ax_bars.spines['top'].set_visible(False)
    ax_bars.spines['right'].set_visible(False)
    ax_bars.grid(axis='y', alpha=0.3, linestyle='--')
    
    # Set sizes bar chart (bottom left)
    ax_sets = fig.add_subplot(gs[2, 0])
    set_sizes = [len(id_dict[key]) for key in tool_keys]
    y_pos = np.arange(n_tools)
    ax_sets.barh(y_pos, set_sizes, color='#555555', edgecolor='black', linewidth=0.8)
    
    # Add set size labels
    for i, size in enumerate(set_sizes):
        ax_sets.text(size + max(set_sizes)*0.02, i, f'{size}', 
                    ha='left', va='center', fontsize=11, fontweight='bold')
    
    ax_sets.set_yticks(y_pos)
    ax_sets.set_yticklabels(tools, fontsize=12, fontweight='bold')
    ax_sets.set_xlabel('Set Size', fontsize=12, fontweight='bold')
    ax_sets.set_xlim(0, max(set_sizes) * 1.15)
    ax_sets.invert_yaxis()
    ax_sets.spines['top'].set_visible(False)
    ax_sets.spines['right'].set_visible(False)
    ax_sets.grid(axis='x', alpha=0.3, linestyle='--')
    
    # Intersection matrix (bottom right)
    ax_matrix = fig.add_subplot(gs[2, 1])
    
    # Draw connection dots and lines
    for combo_idx, combo in enumerate(combinations):
        active_tools = [i for i, val in enumerate(combo) if val == 1]
        
        # Draw dots
        for tool_idx in range(n_tools):
            if tool_idx in active_tools:
                circle = patches.Circle((combo_idx, tool_idx), 0.15, 
                                       facecolor='#313695', edgecolor='black', linewidth=1.5, zorder=3)
            else:
                circle = patches.Circle((combo_idx, tool_idx), 0.1, 
                                       facecolor='lightgray', edgecolor='gray', linewidth=0.5, zorder=2)
            ax_matrix.add_patch(circle)
        
        # Draw connecting line
        if len(active_tools) > 1:
            y_coords = active_tools
            ax_matrix.plot([combo_idx]*len(y_coords), y_coords, 
                          color='#313695', linewidth=3, zorder=1)
    
    ax_matrix.set_xlim(-0.5, n_combos - 0.5)
    ax_matrix.set_ylim(-0.5, n_tools - 0.5)
    ax_matrix.set_yticks(range(n_tools))
    ax_matrix.set_yticklabels([])
    ax_matrix.set_xticks([])
    ax_matrix.invert_yaxis()
    ax_matrix.set_aspect('equal')
    ax_matrix.axis('off')
    
    # Title
    fig.suptitle(title, fontsize=18, fontweight='bold', y=0.98)
    
    # Add note about GT percentage
    fig.text(0.5, 0.02, f'Percentages show proportion of Ground Truth (n={gt_total})', 
            ha='center', fontsize=11, style='italic', color='gray')
    
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"✓ Saved: {output_file}")

# Generate both plots
create_upset_plot(tp_ids, 'Tool Overlap: True Positives Only (TP)', 'combined_upset_plot.png', gt_total)

print("✓ UpSet analysis complete")
PYUPSET
  
  cd - > /dev/null
  echo "  ✓ UpSet plot complete"
fi

################################################################################
# FINAL SUMMARY
################################################################################
log_step "Pipeline Complete"

echo ""
echo "Output structure:"
echo "  Data:  ${OUTPUT_DIR}/data/"
echo "  Plots: ${OUTPUT_DIR}/plots/"
echo ""
echo "Generated files:"
find "${OUTPUT_DIR}/plots" -name "*.png" -type f | wc -l | xargs echo "  PNG plots:"
find "${OUTPUT_DIR}/data" -name "*.csv" -type f | wc -l | xargs echo "  CSV files:"
echo ""
echo "Key plots generated:"
echo "  Length:         ${OUTPUT_DIR}/plots/length/length_barplot.png"
echo "  Classification: ${OUTPUT_DIR}/plots/classification/"
echo "  Performance:    ${OUTPUT_DIR}/plots/performance/tool_performance.png"
echo "                  ${OUTPUT_DIR}/plots/performance/tool_contingency_heatmap.png"
echo "  Combos:         ${OUTPUT_DIR}/plots/combos/combinations_performance.png"
echo "  Expression:     ${OUTPUT_DIR}/plots/expression/expression_correlation.png"
echo "  UpSet:          ${OUTPUT_DIR}/plots/upset/combined_upset_plot.png"
echo ""
echo "=============================================================================="