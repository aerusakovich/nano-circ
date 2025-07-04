#!/usr/bin/env python3
import os
import pandas as pd
from collections import defaultdict
import subprocess
import logging
import matplotlib.pyplot as plt

# Default paths
GTF_FILE = "/home/genouest/cnrs_umr6290/aerusakovich/benchmark/data/ciri-long/files/m38/gencode.vM10.sorted.annotation.gtf"
BEDTOOLS_CMD = "bedtools"

TOOL_BEDS = {
    "ground_truth": "/scratch/aerusakovich/sim_ciri_long_jobim/combined/circRNAs.bed",
    "ciri_long": "/scratch/aerusakovich/sim_ciri_long_jobim/all/bed12/ciri-long/ciri_long_output/CIRI-long.cand_circ.bed12",
    "isocirc": "/scratch/aerusakovich/sim_ciri_long_jobim/combined/isocirc_output/isocirc.bed",
    "circnick": "/scratch/aerusakovich/sim_ciri_long_jobim/all/bed12/circnick/combined/comprehensive_circrna.bed12"
}

TEMP_DIR = "temp_classification"

CIRCRNA_COLORS = {
    'eciRNA': '#1b9e77',    # Teal
    'EIciRNA': '#d95f02',   # Orange
    'ciRNA': '#7570b3',     # Purple-blue
    'intergenic': '#e7298a' # Magenta
}

def run(cmd):
    result = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if result.returncode != 0:
        logging.warning(f"Command failed: {cmd}\n{result.stderr}")
    return result.stdout

def prepare_feature_beds(gtf_file, temp_dir):
    os.makedirs(temp_dir, exist_ok=True)
    gene_bed = os.path.join(temp_dir, "genes.bed")
    exon_bed = os.path.join(temp_dir, "exons.bed")

    if not os.path.exists(gene_bed):
        os.system(f"awk '$3 == \"gene\"' {gtf_file} | awk '{{OFS=\"\\t\"; print $1, $4-1, $5, \".\", 0, $7}}' > {gene_bed}")
    if not os.path.exists(exon_bed):
        os.system(f"awk '$3 == \"exon\"' {gtf_file} | awk '{{OFS=\"\\t\"; print $1, $4-1, $5, \".\", 0, $7}}' > {exon_bed}")

    return gene_bed, exon_bed

def parse_ground_truth_labels(gt_path):
    label_map = {}
    df = pd.read_csv(gt_path, sep='\t', header=None, comment='#')
    if df.shape[1] < 6:
        for i in range(df.shape[1], 6):
            df[i] = '+' if i == 5 else '.'
    df = df[[0, 1, 2, 3]]
    df.columns = ['chrom', 'start', 'end', 'name']

    for name in df['name']:
        parts = name.split('|')
        if len(parts) > 2:
            label_map[name] = parts[2]
    return label_map

def classify_batch(bed_path, bedtools_cmd, gene_bed, exon_bed, temp_dir, ground_truth_labels=None, tool_name="tool"):
    df = pd.read_csv(bed_path, sep='\t', header=None, comment='#')
    if df.shape[1] < 6:
        for i in range(df.shape[1], 6):
            df[i] = '+' if i == 5 else '.'
    df = df[[0, 1, 2, 3, 5]]
    df.columns = ['chrom', 'start', 'end', 'name', 'strand']

    temp_bed = os.path.join(temp_dir, f"{tool_name}_combined.bed")
    df.to_csv(temp_bed, sep='\t', header=False, index=False)

    def extract_ids(cmd):
        return set(line.split('\t')[3] for line in run(cmd).strip().split('\n') if line.strip())

    gene_s = extract_ids(f"{bedtools_cmd} intersect -a {temp_bed} -b {gene_bed} -s -wa -wb")
    exon_s = extract_ids(f"{bedtools_cmd} intersect -a {temp_bed} -b {exon_bed} -s -wa -wb")
    full_exon_s = extract_ids(f"{bedtools_cmd} intersect -a {temp_bed} -b {exon_bed} -s -f 1.0 -wa")

    gene_u = extract_ids(f"{bedtools_cmd} intersect -a {temp_bed} -b {gene_bed} -wa -wb")
    exon_u = extract_ids(f"{bedtools_cmd} intersect -a {temp_bed} -b {exon_bed} -wa -wb")
    full_exon_u = extract_ids(f"{bedtools_cmd} intersect -a {temp_bed} -b {exon_bed} -f 1.0 -wa")

    counts = defaultdict(int)
    intergenic_coords = []

    for _, row in df.iterrows():
        name = row['name']
        chrom, start, end, strand = row['chrom'], row['start'], row['end'], row['strand']

        in_gene = name in gene_s or name in gene_u
        in_full_exon = name in full_exon_s or name in full_exon_u
        in_exon = name in exon_s or name in exon_u

        if in_gene:
            if in_full_exon:
                circ_type = "eciRNA"
            elif in_exon:
                circ_type = "EIciRNA"
            else:
                circ_type = "ciRNA"
        else:
            circ_type = "intergenic"
            gt_label = ground_truth_labels.get(name, "not_in_gt")
            print(f"[DEBUG: {tool_name}] {name} classified as intergenic (GT label: {gt_label})")
            intergenic_coords.append((chrom, start, end, name, '.', strand))

        counts[circ_type] += 1

    if intergenic_coords:
        intergenic_bed = os.path.join(temp_dir, f"{tool_name}_intergenic_coords.bed")
        with open(intergenic_bed, 'w') as f:
            for coord in intergenic_coords:
                f.write('\t'.join(map(str, coord)) + '\n')
        print(f"Saved intergenic coordinates to: {intergenic_bed}")

    os.remove(temp_bed)
    return counts

def create_improved_pie_chart(data, labels, colors, title, filename):
    total = sum(data)
    percentages = [100 * val / total for val in data]

    plt.figure(figsize=(16, 10))
    gs = plt.GridSpec(1, 2, width_ratios=[2, 1])
    ax_pie = plt.subplot(gs[0])
    wedges, _ = ax_pie.pie(
        data,
        colors=[colors[label] for label in labels],
        startangle=90,
        wedgeprops={'linewidth': 1.5, 'edgecolor': 'black'},
        radius=0.9
    )

    legend_labels = [
        f"{label}: {percent:.1f}% ({val})"
        for val, label, percent in zip(data, labels, percentages)
    ]

    ax_legend = plt.subplot(gs[1])
    ax_legend.axis('off')
    legend_elements = [
        plt.Rectangle((0, 0), 1, 1, color=colors[labels[i]], edgecolor='black', linewidth=1.5)
        for i in range(len(labels))
    ]
    ax_legend.legend(
        legend_elements,
        legend_labels,
        loc='center left',
        title=title,
        fontsize=16,
        title_fontsize=20,
        frameon=True,
        edgecolor='black',
        facecolor='white',
        framealpha=1,
        bbox_to_anchor=(0, 0.5)
    )

    plt.tight_layout()
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close()

def plot_results(all_results):
    out_dir = os.path.join(TEMP_DIR, "plots")
    os.makedirs(out_dir, exist_ok=True)

    for tool, result in all_results.items():
        labels = [k for k in ["eciRNA", "EIciRNA", "ciRNA", "intergenic"] if k in result]
        data = [result[k] for k in labels]
        output_path = os.path.join(out_dir, f"{tool}_classification_pie.png")
        create_improved_pie_chart(data, labels, CIRCRNA_COLORS, f"{tool} classification", output_path)

def classify_all_tools(tool_beds, gtf_file, bedtools_cmd, temp_dir):
    gene_bed, exon_bed = prepare_feature_beds(gtf_file, temp_dir)
    ground_truth_labels = parse_ground_truth_labels(tool_beds["ground_truth"])

    all_results = {}

    for tool_name, bed_path in tool_beds.items():
        print(f"\nClassifying {tool_name} from {bed_path}")
        if not os.path.exists(bed_path):
            print(f"  ❌ Skipping (file not found)")
            continue

        counts = classify_batch(
            bed_path,
            bedtools_cmd,
            gene_bed,
            exon_bed,
            temp_dir,
            ground_truth_labels=ground_truth_labels,
            tool_name=tool_name
        )

        print(f"  ✅ Classified {sum(counts.values())} circRNAs:")
        for t in ['eciRNA', 'EIciRNA', 'ciRNA', 'intergenic']:
            if t in counts:
                print(f"    {t}: {counts[t]}")
        all_results[tool_name] = dict(counts)

    plot_results(all_results)
    print(f"\nAll results and plots saved to: {os.path.abspath(temp_dir)}")

    return all_results

if __name__ == "__main__":
    classify_all_tools(TOOL_BEDS, GTF_FILE, BEDTOOLS_CMD, TEMP_DIR)
