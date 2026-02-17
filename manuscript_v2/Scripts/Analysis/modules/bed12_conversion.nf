process CONVERT_TO_BED12_CIRILONG {
    tag "${run_name}"
    publishDir "${params.outdir}/bed12/cirilong/${run_name}", mode: 'copy'
    
    input:
    tuple val(run_name), path(output_dir), path(metrics)
    
    output:
    tuple val(run_name), path("CIRI-long.cand_circ.bed12"), emit: bed12
    
    script:
    """
    #!/usr/bin/env python3
    
    import os
    import re
    try:
        from Bio import SeqIO
    except ImportError:
        import subprocess
        subprocess.check_call(["pip", "install", "biopython", "--user", "--quiet"])
        from Bio import SeqIO
    
    def parse_ciri_info(info_file):
        circ_info = {}
        info_count = 0
        multi_exon_count = 0
        
        with open(info_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                fields = line.strip().split('\\t')
                if len(fields) < 9:
                    continue
                    
                info_count += 1
                chrom = fields[0]
                start = int(fields[3])  # 1-based in GFF
                end = int(fields[4])
                score = fields[5]
                strand = fields[6]
                
                # Extract attributes
                attr_dict = {}
                attrs = fields[8].split(';')
                for attr in attrs:
                    attr = attr.strip()
                    if not attr:
                        continue
                    if '"' in attr:
                        key, value = attr.split(' "', 1)
                        value = value.strip('"')
                        attr_dict[key.strip()] = value
                    else:
                        try:
                            key, value = attr.split(' ', 1)
                            attr_dict[key.strip()] = value.strip()
                        except ValueError:
                            continue
                
                circ_id = attr_dict.get('circ_id', '').strip('"')
                if not circ_id:
                    circ_id = f"{chrom}:{start}-{end}"
                
                # Parse isoform information
                isoform = attr_dict.get('isoform', '').strip('"')
                exons = []
                
                if isoform:
                    exon_coords = isoform.split(',')
                    for coord in exon_coords:
                        if '-' in coord:
                            try:
                                parts = coord.split('-')
                                if len(parts) == 2:
                                    exon_start = int(parts[0])
                                    exon_end = int(parts[1])
                                    exons.append((exon_start, exon_end))
                            except ValueError as e:
                                continue
                
                # If no exon information, treat as single exon
                if not exons:
                    exons = [(start, end)]
                
                if len(exons) > 1:
                    multi_exon_count += 1
                    
                circ_info[circ_id] = {
                    'chrom': chrom,
                    'start': start,
                    'end': end,
                    'score': score,
                    'strand': strand,
                    'exons': sorted(exons, key=lambda x: x[0]),
                    'circ_id': circ_id,
                    'gene_name': attr_dict.get('gene_name', '').strip('"'),
                    'isoform': isoform
                }
        
        print(f"Parsed {info_count} entries from CIRI-long info file, extracted {len(circ_info)} unique circRNAs")
        print(f"Found {multi_exon_count} multi-exon circRNAs in info file")
        return circ_info
    
    # Find the info file
    info_file = None
    for root, dirs, files in os.walk('${output_dir}'):
        for f in files:
            if f.endswith('.info'):
                info_file = os.path.join(root, f)
                break
        if info_file:
            break
    
    if not info_file:
        print("ERROR: No .info file found")
        exit(1)
    
    print(f"Using info file: {info_file}")
    
    # Parse and convert
    circ_info = parse_ciri_info(info_file)
    
    with open('CIRI-long.cand_circ.bed12', 'w') as out:
        for circ_id, info in circ_info.items():
            chrom = info['chrom']
            start = info['start'] - 1  # Convert to 0-based for BED
            end = info['end']
            name = circ_id
            score = info['score']
            strand = info['strand']
            
            # BED12 requires these fields
            thick_start = start
            thick_end = end
            item_rgb = "0,0,0"
            
            exons = info['exons']
            block_count = len(exons)
            
            # Calculate block sizes and starts (relative to start)
            block_sizes = []
            block_starts = []
            
            for exon_start, exon_end in exons:
                block_size = exon_end - exon_start + 1
                block_start = exon_start - 1 - start
                
                block_sizes.append(str(block_size))
                block_starts.append(str(block_start))
            
            # Join with commas
            block_sizes_str = ','.join(block_sizes)
            block_starts_str = ','.join(block_starts)
            
            # Write BED12 line
            bed_line = f"{chrom}\\t{start}\\t{end}\\t{name}\\t{score}\\t{strand}\\t{thick_start}\\t{thick_end}\\t{item_rgb}\\t{block_count}\\t{block_sizes_str}\\t{block_starts_str}\\n"
            out.write(bed_line)
    
    print(f"Successfully created BED12 file with {len(circ_info)} entries")
    """
}

process CONVERT_TO_BED12_CIRCNICK {
    tag "${run_name}"
    publishDir "${params.outdir}/bed12/circnick/${run_name}", mode: 'copy'
    
    input:
    tuple val(run_name), path(output_dir), path(metrics)
    
    output:
    tuple val(run_name), path("comprehensive_circrna.bed12"), emit: bed12
    
    script:
    """
    #!/usr/bin/env python3
    
    import os
    import sys
    
    def parse_exons(exons_file):
        exons_by_chrom = {}
        with open(exons_file, 'r') as f:
            for line in f:
                try:
                    parts = line.strip().split('\\t')
                    chrom = parts[0]
                    start = int(parts[1])
                    end = int(parts[2])
                    strand = parts[5] if len(parts) > 5 else '+'
                    
                    if chrom not in exons_by_chrom:
                        exons_by_chrom[chrom] = []
                    exons_by_chrom[chrom].append((start, end, strand))
                except (ValueError, IndexError):
                    continue
        
        for chrom in exons_by_chrom:
            exons_by_chrom[chrom].sort(key=lambda x: x[0])
        
        return exons_by_chrom
    
    def find_matching_exons(chrom, start, end, strand, exons_by_chrom):
        matching_exons = []
        
        if chrom not in exons_by_chrom:
            return [(start, end)]
        
        for ex_start, ex_end, ex_strand in exons_by_chrom[chrom]:
            if ex_strand != strand:
                continue
                
            if (ex_start >= start and ex_end <= end) or \\
               (ex_start < start and ex_end > start) or \\
               (ex_start < end and ex_end > end):
                effective_start = max(start, ex_start)
                effective_end = min(end, ex_end)
                matching_exons.append((effective_start, effective_end))
        
        if not matching_exons:
            matching_exons.append((start, end))
        
        matching_exons.sort(key=lambda x: x[0])
        merged_exons = []
        
        if matching_exons:
            current = matching_exons[0]
            for next_ex in matching_exons[1:]:
                if current[1] >= next_ex[0]:
                    current = (current[0], max(current[1], next_ex[1]))
                else:
                    merged_exons.append(current)
                    current = next_ex
            merged_exons.append(current)
        
        return merged_exons
    
    # Find candidates and exons files
    candidates_file = None
    all_exons_file = None
    
    for root, dirs, files in os.walk('${output_dir}'):
        for f in files:
            if 'circRNA_candidates.annotated.txt' in f:
                candidates_file = os.path.join(root, f)
            if 'allExons.bed' in f and 'notGencode' not in f:
                all_exons_file = os.path.join(root, f)
    
    if not candidates_file or not all_exons_file:
        print("ERROR: Required files not found")
        exit(1)
    
    print(f"Using candidates: {candidates_file}")
    print(f"Using exons: {all_exons_file}")
    
    # Parse exons
    exons_by_chrom = parse_exons(all_exons_file)
    
    # Process candidates
    with open('comprehensive_circrna.bed12', 'w') as out:
        with open(candidates_file, 'r') as f:
            next(f)  # Skip header
            
            for line in f:
                try:
                    parts = line.strip().split('\\t')
                    if len(parts) < 7:
                        continue
                        
                    circ_id = parts[0]
                    chrom = parts[1]
                    start = int(parts[2])
                    end = int(parts[3])
                    strand = parts[6]
                    
                    if start >= end or strand not in ['+', '-']:
                        continue
                        
                except (ValueError, IndexError):
                    continue
                
                # Find matching exons
                matching_exons = find_matching_exons(chrom, start, end, strand, exons_by_chrom)
                
                exon_sizes = []
                exon_starts = []
                
                for ex_start, ex_end in matching_exons:
                    exon_sizes.append(str(ex_end - ex_start))
                    exon_starts.append(str(ex_start - start))
                
                if not exon_sizes:
                    continue
                
                bed_entry = [
                    chrom, 
                    str(start), 
                    str(end), 
                    circ_id, 
                    '1000', 
                    strand, 
                    str(start), 
                    str(end), 
                    '0,0,0', 
                    str(len(matching_exons)), 
                    ','.join(exon_sizes) + ',',
                    ','.join(exon_starts) + ','
                ]
                
                if len(bed_entry) == 12:
                    out.write('\\t'.join(bed_entry) + '\\n')
    
    print("CircNick BED12 conversion complete")
    """
}

process CONVERT_TO_BED12_GROUNDTRUTH {
    tag "${run_name}"
    publishDir "${params.outdir}/bed12/groundtruth/${run_name}", mode: 'copy'

    cpus = params.threads
    memory = '16 GB'
    time   = '4h'
    clusterOptions = ''

    input:
    tuple val(run_name), path(run_dir)
    path genome

    output:
    tuple val(run_name), path("ground_truth_circular.bed12"), emit: bed12_circ
    tuple val(run_name), path("ground_truth.sorted.bam"),     emit: bam
    tuple val(run_name), path("ground_truth.sorted.bam.bai"), emit: bai
    path "filtering_summary.txt"
    path "ground_truth_counts.debug.tsv"

    script:
    // Groovy-side defaults (override in params.yaml)
    def min_reads             = params.min_reads_per_circ    ?: 1
    def mapq_threshold        = params.mapq_threshold        ?: 0
    def initial_tile          = params.initial_tile          ?: 1
    def rescue_tile           = params.rescue_tile           ?: 5
    def rescue_enabled        = params.rescue_enabled        ?: 1
    def rescue_opts           = params.rescue_opts           ?: '-k11 -w2 -N50'
    def mode                  = params.mode                  ?: 'best'    // best | unique_any
    def rescue_preset         = params.rescue_preset         ?: ''       // "" = no preset
    def rescue_mapq           = params.rescue_mapq           ?: 0
    def rescue_expected_ratio = params.rescue_expected_ratio ?: 0.25     // alpha
    def sim_dir               = params.simulated_data_dir    ?: ''

    // BLAT-specific params
    def blat_enabled          = params.blat_enabled          ?: 0
    def blat_tile             = params.blat_tile             ?: 10
    // IMPORTANT: no -fastMap by default (it breaks >5kb)
    def blat_opts             = params.blat_opts             ?: '-minScore=20 -minIdentity=70'

    """
    #!/bin/bash

    source /local/env/envpython-3.9.5.sh || true
    conda activate bed12 || true

    set -euo pipefail

    echo "Run: ${run_name}"
    echo "Genome: ${genome}"
    echo "Simulation dir: ${sim_dir}"

    ########################################
    # 0) Locate reads and original BED
    ########################################

    READS=""
    if [ -f "${run_dir}/pooled/combined_reads.fastq" ]; then
        READS="${run_dir}/pooled/combined_reads.fastq"
    elif [ -f "${run_dir}/pooled/combined_reads.fq.gz" ]; then
        READS="${run_dir}/pooled/combined_reads.fq.gz"
    fi

    if [ -z "\$READS" ]; then
        echo "ERROR: No combined_reads found in ${run_dir}/pooled" >&2
        exit 1
    fi

    if [ -z "${sim_dir}" ]; then
        echo "ERROR: params.simulated_data_dir not set" >&2
        exit 1
    fi

    ORIG_BED_RAW="${sim_dir}/${run_name}/circRNAs.bed"
    if [ ! -f "\$ORIG_BED_RAW" ]; then
        echo "ERROR: Original circRNAs.bed not found at \$ORIG_BED_RAW" >&2
        exit 1
    fi

    ABUND_FILE="${sim_dir}/${run_name}/abundances.tsv"
    export ABUND_FILE

    HAVE_ABUND=0
    if [ -f "\$ABUND_FILE" ]; then
        HAVE_ABUND=1
    else
        echo "WARNING: abundances.tsv not found at \$ABUND_FILE – TPM-based rescue will be skipped."
    fi

    # Deduplicate BED12 by name (keep first occurrence)
    awk '!seen[\$4]++' "\$ORIG_BED_RAW" > circRNAs_dedup.bed
    export ORIG_BED="circRNAs_dedup.bed"

    TOTAL_ORIG=\$(wc -l < "\$ORIG_BED_RAW")
    TOTAL_DEDUP=\$(wc -l < "\$ORIG_BED")
    DUPLICATES=\$((TOTAL_ORIG - TOTAL_DEDUP))

    echo "Original BED12:  \$TOTAL_ORIG entries"
    echo "After dedup:     \$TOTAL_DEDUP entries (\$DUPLICATES duplicates removed)"
    echo "Using reads:     \$READS"
    echo "Using ORIG_BED:  \$ORIG_BED"
    echo "Parameters:"
    echo "  MAPQ threshold: ${mapq_threshold}"
    echo "  Mode: ${mode}"
    echo "  Initial tiling: ${initial_tile}x"
    echo "  Rescue tiling:  ${rescue_tile}x"
    echo "  Rescue enabled: ${rescue_enabled}"
    echo "  Rescue opts:    ${rescue_opts}"
    echo "  Rescue preset:  '${rescue_preset}' (empty = no preset)"
    echo "  Rescue MAPQ:    ${rescue_mapq}"
    echo "  Rescue alpha:   ${rescue_expected_ratio}"
    echo "  Min reads:      ${min_reads}"
    echo "  BLAT enabled:   ${blat_enabled}"
    echo "  BLAT tiling:    ${blat_tile}x"
    echo "  BLAT opts:      ${blat_opts}"

    ########################################
    # 1) Extract circ-only reads
    ########################################

    if [[ "\$READS" =~ \\.gz\$ ]]; then
        zcat "\$READS" | awk 'NR%4==1{h=\$0; keep=(h ~ /_circ/)} keep' > circ_only.fastq
    else
        awk 'NR%4==1{h=\$0; keep=(h ~ /_circ/)} keep' "\$READS" > circ_only.fastq
    fi

    CIRC_N=\$(awk 'END{print NR/4}' circ_only.fastq)
    echo "Found \$CIRC_N circ reads in FASTQ"

    ########################################
    # 2) Header counts by key (from FULL reads)
    ########################################

    if [[ "\$READS" =~ \\.gz\$ ]]; then
        zcat "\$READS" | \\
        awk 'NR%4==1{
                 h=\$0; sub(/^@/,"",h); sub(/ .*/,"",h);
                 key=h; sub(/\\|.*/,"",key); sub(/_.*/,"",key);
                 cnt[key]++
             }
             END{
                 for(k in cnt) print k"\\t"cnt[k]
             }' > header_counts.tsv
    else
        awk 'NR%4==1{
                 h=\$0; sub(/^@/,"",h); sub(/ .*/,"",h);
                 key=h; sub(/\\|.*/,"",key); sub(/_.*/,"",key);
                 cnt[key]++
             }
             END{
                 for(k in cnt) print k"\\t"cnt[k]
             }' "\$READS" > header_counts.tsv
    fi

    ########################################
    # 3) Build circRNA-tiled FASTA from BED12 + genome
    ########################################

    samtools faidx "${genome}" || true

    python3 << 'EOFPYTHON'
import subprocess, os

bed_file = os.environ.get('ORIG_BED')
genome   = "${genome}"
tile_n   = int(${initial_tile})

if not bed_file or not os.path.exists(bed_file):
    raise SystemExit("ORIG_BED not set or missing")

def bed12_to_fasta(bed_file, genome_file, output_fasta, tile_n):
    with open(bed_file) as bed, open(output_fasta, "w") as fa:
        for line in bed:
            if line.startswith("#") or line.startswith("track"):
                continue
            fields = line.rstrip().split()
            if len(fields) < 12:
                continue
            chrom = fields[0]
            start = int(fields[1])
            name  = fields[3]
            strand = fields[5]
            block_count = int(fields[9])
            block_sizes  = [int(x) for x in fields[10].rstrip(",").split(",")]
            block_starts = [int(x) for x in fields[11].rstrip(",").split(",")]

            exon_seqs = []
            for i in range(block_count):
                exon_start = start + block_starts[i]
                exon_end   = exon_start + block_sizes[i]
                cmd = f"samtools faidx {genome_file} {chrom}:{exon_start+1}-{exon_end}"
                res = subprocess.run(cmd, shell=True, capture_output=True, text=True, check=True)
                seq = "".join(res.stdout.splitlines()[1:])
                exon_seqs.append(seq)

            full_seq = "".join(exon_seqs)
            if strand == "-":
                comp = str.maketrans("ATCGatcg", "TAGCtagc")
                full_seq = full_seq.translate(comp)[::-1]

            print(f">{name}", file=fa)
            print(full_seq * tile_n, file=fa)

bed12_to_fasta(bed_file, genome, "circRNAs_tiled.fa", tile_n)
print(f"Created circRNAs_tiled.fa with {tile_n}x tiling")
EOFPYTHON

    ########################################
    # 4) Initial mapping & counting (minimap2)
    ########################################

    echo "Mapping circ reads to circRNA sequences (initial, map-ont)..."
    minimap2 -ax map-ont -t ${task.cpus} circRNAs_tiled.fa circ_only.fastq 2> minimap2.initial.log | \\
        samtools view -bS - | samtools sort -@ ${task.cpus} -o initial.bam -
    samtools index initial.bam

    echo "Counting initial mappings (MODE='${mode}', MAPQ >= ${mapq_threshold})..."

    if [ "${mode}" = "best" ]; then
        samtools view -q ${mapq_threshold} -F 2308 initial.bam | \\
            awk '{print \$3}' | sort | uniq -c | \\
            awk '{print \$2"\\t"\$1}' > initial_counts.tsv
    else
        # unique_any: keep only reads mapping to exactly one circRNA at >= MAPQ
        samtools view -q ${mapq_threshold} initial.bam | \\
          awk '
            { q=\$1; r=\$3; if(r=="*"||r=="") next;
              key=q SUBSEP r;
              if(!(key in seen)){ seen[key]=1; hit[q]++; last[r,q]=1 }
            }
            END{
              for (q in hit) if (hit[q]==1) {
                for (rk in last) {
                  split(rk, a, SUBSEP)
                  if (a[2]==q) print a[1]
                }
              }
            }' | \\
          sort | uniq -c | awk '{print \$2"\\t"\$1}' > initial_counts.tsv
    fi

    INITIAL_MATCHED=\$(wc -l < initial_counts.tsv || echo 0)
    echo "Initial mapping: \$INITIAL_MATCHED circRNAs with >=1 read"

    ########################################
    # 5) id -> key map from ORIG_BED
    ########################################

    awk 'BEGIN{FS=OFS="\\t"}
         {
           id=\$4;
           key=id;
           sub(/\\|.*/,"",key);
           sub(/\\..*/,"",key);
           print id, key
         }' "\$ORIG_BED" > id2key.tsv

    ########################################
    # 6) TPM -> expected_reads (if abundances.tsv exists)
    ########################################

    if [ "\$HAVE_ABUND" -eq 1 ]; then
        python3 << 'EOFPYTHON'
import os

ab_file = os.environ.get("ABUND_FILE")
if not ab_file or not os.path.exists(ab_file):
    raise SystemExit("ABUND_FILE missing")

tpm = {}
total_tpm = 0.0

with open(ab_file) as f:
    header = f.readline().rstrip().split()
    hl = [h.lower() for h in header]
    try:
        id_idx = hl.index("circrna_id")
    except ValueError:
        try:
            id_idx = hl.index("id")
        except ValueError:
            id_idx = 0
    try:
        tpm_idx = hl.index("tpm")
    except ValueError:
        raise SystemExit("No TPM column in abundances.tsv")

    for line in f:
        if not line.strip():
            continue
        cols = line.rstrip().split()
        if len(cols) <= max(id_idx, tpm_idx):
            continue
        cid = cols[id_idx]
        vs  = cols[tpm_idx]
        try:
            v = float(vs) if vs else 0.0
        except ValueError:
            v = 0.0
        tpm[cid] = v
        total_tpm += v

if total_tpm <= 0.0:
    total_tpm = 1.0

# count circ reads
with open("circ_only.fastq") as f:
    lines = sum(1 for _ in f)
circ_n = lines // 4 if lines else 0

with open("expected_reads.tsv", "w") as out:
    print("id\texpected_reads", file=out)
    for cid, v in tpm.items():
        exp = v * circ_n / total_tpm
        print(f"{cid}\t{exp:.3f}", file=out)
EOFPYTHON
    else
        echo -e "id\\texpected_reads" > expected_reads.tsv
    fi

    ########################################
    # 7) Define rescue candidates (minimap2 rescue)
    ########################################

    if [ ${rescue_enabled} -eq 1 ]; then
        echo "Defining rescue candidates (header + TPM-aware)..."
        python3 << 'EOFPYTHON'
import os
from collections import defaultdict

alpha = float(${rescue_expected_ratio})

# initial counts
mapped = {}
if os.path.exists("initial_counts.tsv"):
    with open("initial_counts.tsv") as f:
        for line in f:
            if not line.strip():
                continue
            cid, cnt = line.rstrip().split()
            mapped[cid] = int(cnt)

# header counts by key
headers = defaultdict(int)
if os.path.exists("header_counts.tsv"):
    with open("header_counts.tsv") as f:
        for line in f:
            if not line.strip():
                continue
            k, c = line.rstrip().split()
            headers[k] = int(c)

# expected reads by id
expected = {}
if os.path.exists("expected_reads.tsv"):
    with open("expected_reads.tsv") as f:
        next(f, None)
        for line in f:
            if not line.strip():
                continue
            cid, e = line.rstrip().split()
            try:
                expected[cid] = float(e)
            except ValueError:
                expected[cid] = 0.0

def id_to_key(cid: str) -> str:
    k = cid.split("|", 1)[0]
    k = k.split(".", 1)[0]
    return k

rescue = set()

orig_bed = os.environ.get("ORIG_BED")
if not orig_bed or not os.path.exists(orig_bed):
    raise SystemExit("ORIG_BED missing")

with open(orig_bed) as f:
    for line in f:
        if line.startswith("#") or line.startswith("track"):
            continue
        fields = line.rstrip().split()
        if len(fields) < 4:
            continue
        cid = fields[3]
        key = id_to_key(cid)
        obs = mapped.get(cid, 0)
        hdr = headers.get(key, 0)
        exp = expected.get(cid, 0.0)

        # 1) header>0 but obs==0
        if hdr > 0 and obs == 0:
            rescue.add(cid)
            continue
        # 2) exp>0 but obs==0
        if exp > 0 and obs == 0:
            rescue.add(cid)
            continue
        # 3) obs << expected (for reasonably high expected)
        if exp >= 2.0 and obs < alpha * exp:
            rescue.add(cid)

with open("rescue_ids.txt", "w") as out:
    for cid in sorted(rescue):
        print(cid, file=out)

print(f"Rescue candidates: {len(rescue)}")
EOFPYTHON
    else
        echo "Rescue disabled; creating empty rescue_ids.txt"
        > rescue_ids.txt
    fi

    RESCUE_N=\$(wc -l < rescue_ids.txt || echo 0)

    ########################################
    # 8) Extract rescue reads via key-table (minimap2 rescue)
    ########################################

    RESCUE_READS=0
    if [ ${rescue_enabled} -eq 1 ] && [ "\$RESCUE_N" -gt 0 ]; then
        echo "Preparing rescue reads for \$RESCUE_N circRNAs..."

        awk 'FNR==NR{ id2key[\$1]=\$2; next }
             (\$1 in id2key){ print id2key[\$1] }' id2key.tsv rescue_ids.txt | sort -u > rescue.keys

        awk -v K="rescue.keys" '
          BEGIN{
            while((getline line < K) > 0){ keys[line]=1 }
          }
          NR%4==1{
            h=\$0; sub(/^@/,"",h); sub(/ .*/,"",h);
            key=h; sub(/\\|.*/,"",key); sub(/_.*/,"",key);
            keep=(key in keys)
          }
          { if(keep) print }
        ' circ_only.fastq > rescue.fastq

        RESCUE_READS=\$(awk 'END{print NR/4}' rescue.fastq || echo 0)
        echo "Extracted \$RESCUE_READS rescue reads"
    else
        echo "No rescue candidates or rescue disabled"
        > rescue.keys
        > rescue.fastq
    fi

    ########################################
    # 9) Build rescue FASTA & run rescue mapping (minimap2)
    ########################################

    > rescue_counts.tsv
    if [ ${rescue_enabled} -eq 1 ] && [ "\$RESCUE_N" -gt 0 ] && [ "\$RESCUE_READS" -gt 0 ]; then
        echo "Building rescue FASTA with ${rescue_tile}x tiling..."
        python3 << 'EOFPYTHON'
import subprocess, os

bed_file = os.environ.get('ORIG_BED')
genome   = "${genome}"
tile_n   = int(${rescue_tile})

if not bed_file or not os.path.exists(bed_file):
    raise SystemExit("ORIG_BED missing")

rescue_ids = set()
with open("rescue_ids.txt") as f:
    for line in f:
        line = line.strip()
        if line:
            rescue_ids.add(line)

def bed12_to_fasta_rescue(bed_file, genome_file, output_fasta, tile_n, rescue_ids):
    with open(bed_file) as bed, open(output_fasta, "w") as fa:
        for line in bed:
            if line.startswith("#") or line.startswith("track"):
                continue
            fields = line.rstrip().split()
            if len(fields) < 12:
                continue
            name = fields[3]
            if name not in rescue_ids:
                continue
            chrom = fields[0]
            start = int(fields[1])
            strand = fields[5]
            block_count = int(fields[9])
            block_sizes  = [int(x) for x in fields[10].rstrip(",").split(",")]
            block_starts = [int(x) for x in fields[11].rstrip(",").split(",")]

            exon_seqs = []
            for i in range(block_count):
                exon_start = start + block_starts[i]
                exon_end   = exon_start + block_sizes[i]
                cmd = f"samtools faidx {genome_file} {chrom}:{exon_start+1}-{exon_end}"
                res = subprocess.run(cmd, shell=True, capture_output=True, text=True, check=True)
                seq = "".join(res.stdout.splitlines()[1:])
                exon_seqs.append(seq)

            full_seq = "".join(exon_seqs)
            if strand == "-":
                comp = str.maketrans("ATCGatcg", "TAGCtagc")
                full_seq = full_seq.translate(comp)[::-1]

            print(f">{name}", file=fa)
            print(full_seq * tile_n, file=fa)

bed12_to_fasta_rescue(bed_file, genome, "circRNAs_rescue.fa", tile_n, rescue_ids)
print(f"Created circRNAs_rescue.fa for {len(rescue_ids)} circRNAs (tile {tile_n}x)")
EOFPYTHON

        echo "Running rescue mapping..."
        if [ -s circRNAs_rescue.fa ] && [ -s rescue.fastq ]; then
            if [ -n "${rescue_preset}" ]; then
                minimap2 -ax ${rescue_preset} ${rescue_opts} -t ${task.cpus} circRNAs_rescue.fa rescue.fastq 2> minimap2.rescue.log | \\
                    samtools view -bS - | samtools sort -@ ${task.cpus} -o rescue.bam -
            else
                minimap2 -a ${rescue_opts} -t ${task.cpus} circRNAs_rescue.fa rescue.fastq 2> minimap2.rescue.log | \\
                    samtools view -bS - | samtools sort -@ ${task.cpus} -o rescue.bam -
            fi
            samtools index rescue.bam
            samtools view -q ${rescue_mapq} -F 2308 rescue.bam | \\
                awk '{print \$3}' | sort | uniq -c | \\
                awk '{print \$2"\\t"\$1}' > rescue_counts.tsv || > rescue_counts.tsv
        else
            echo "Rescue FASTA or rescue.fastq empty; creating empty rescue.bam"
            samtools view -b -o rescue.bam /dev/null 2>/dev/null || true
            > rescue_counts.tsv
        fi
    else
        echo "Skipping rescue mapping"
        samtools view -b -o rescue.bam /dev/null 2>/dev/null || true
        > rescue_counts.tsv
    fi

    ########################################
    # 9b) BLAT rescue (third pass)
    ########################################

    > blat_counts.tsv

    if [ ${blat_enabled} -eq 1 ]; then
        echo "Preparing BLAT rescue candidates..."

        python3 << 'EOFPYTHON'
import os

# load initial + rescue counts
initial = {}
if os.path.exists("initial_counts.tsv"):
    with open("initial_counts.tsv") as f:
        for line in f:
            if not line.strip():
                continue
            cid, cnt = line.rstrip().split()
            initial[cid] = int(cnt)

rescue = {}
if os.path.exists("rescue_counts.tsv"):
    with open("rescue_counts.tsv") as f:
        for line in f:
            if not line.strip():
                continue
            cid, cnt = line.rstrip().split()
            rescue[cid] = int(cnt)

# header counts by key
headers = {}
if os.path.exists("header_counts.tsv"):
    with open("header_counts.tsv") as f:
        for line in f:
            if not line.strip():
                continue
            key, cnt = line.rstrip().split()
            headers[key] = int(cnt)

# expected reads
expected = {}
if os.path.exists("expected_reads.tsv"):
    with open("expected_reads.tsv") as f:
        next(f, None)
        for line in f:
            if not line.strip():
                continue
            cid, e = line.rstrip().split()
            try:
                expected[cid] = float(e)
            except ValueError:
                expected[cid] = 0.0

def id_to_key(cid: str) -> str:
    k = cid.split("|", 1)[0]
    k = k.split(".", 1)[0]
    return k

orig_bed = os.environ.get("ORIG_BED")
if not orig_bed or not os.path.exists(orig_bed):
    raise SystemExit("ORIG_BED missing")

blat_ids = set()

with open(orig_bed) as f:
    for line in f:
        if line.startswith("#") or line.startswith("track"):
            continue
        fields = line.rstrip().split()
        if len(fields) < 4:
            continue
        cid = fields[3]
        key = id_to_key(cid)
        init_cnt = initial.get(cid, 0)
        resc_cnt = rescue.get(cid, 0)
        total_cnt = init_cnt + resc_cnt
        hdr = headers.get(key, 0)
        exp = expected.get(cid, 0.0)

        # still unmapped after rescue, but evidence exists
        if total_cnt == 0 and (hdr > 0 or exp > 0.0):
            blat_ids.add(cid)

with open("blat_ids.txt", "w") as out:
    for cid in sorted(blat_ids):
        print(cid, file=out)

print(f"BLAT candidates: {len(blat_ids)}")
EOFPYTHON

        BLAT_N=\$(wc -l < blat_ids.txt || echo 0)

        if [ "\$BLAT_N" -gt 0 ]; then
            echo "Building BLAT FASTA with ${blat_tile}x tiling..."

            python3 << 'EOFPYTHON'
import subprocess, os

bed_file = os.environ.get('ORIG_BED')
genome   = "${genome}"
tile_n   = int(${blat_tile})

if not bed_file or not os.path.exists(bed_file):
    raise SystemExit("ORIG_BED missing")

blat_ids = set()
with open("blat_ids.txt") as f:
    for line in f:
        line = line.strip()
        if line:
            blat_ids.add(line)

def bed12_to_fasta_blat(bed_file, genome_file, output_fasta, tile_n, blat_ids):
    with open(bed_file) as bed, open(output_fasta, "w") as fa:
        for line in bed:
            if line.startswith("#") or line.startswith("track"):
                continue
            fields = line.rstrip().split()
            if len(fields) < 12:
                continue
            name = fields[3]
            if name not in blat_ids:
                continue
            chrom = fields[0]
            start = int(fields[1])
            strand = fields[5]
            block_count = int(fields[9])
            block_sizes  = [int(x) for x in fields[10].rstrip(",").split(",")]
            block_starts = [int(x) for x in fields[11].rstrip(",").split(",")]

            exon_seqs = []
            for i in range(block_count):
                exon_start = start + block_starts[i]
                exon_end   = exon_start + block_sizes[i]
                cmd = f"samtools faidx {genome_file} {chrom}:{exon_start+1}-{exon_end}"
                res = subprocess.run(cmd, shell=True, capture_output=True, text=True, check=True)
                seq = "".join(res.stdout.splitlines()[1:])
                exon_seqs.append(seq)

            full_seq = "".join(exon_seqs)
            if strand == "-":
                comp = str.maketrans("ATCGatcg", "TAGCtagc")
                full_seq = full_seq.translate(comp)[::-1]

            print(f">{name}", file=fa)
            print(full_seq * tile_n, file=fa)

bed12_to_fasta_blat(bed_file, genome, "circRNAs_blat.fa", tile_n, blat_ids)
print(f"Created circRNAs_blat.fa for {len(blat_ids)} circRNAs (tile {tile_n}x)")
EOFPYTHON

            echo "Extracting reads for BLAT..."
            awk 'FNR==NR{ id2key[\$1]=\$2; next }
                 (\$1 in id2key){ print id2key[\$1] }' id2key.tsv blat_ids.txt | sort -u > blat.keys

            awk -v K="blat.keys" '
              BEGIN{
                while((getline line < K) > 0){ keys[line]=1 }
              }
              NR%4==1{
                h=\$0; sub(/^@/,"",h); sub(/ .*/,"",h);
                key=h; sub(/\\|.*/,"",key); sub(/_.*/,"",key);
                keep=(key in keys)
              }
              { if(keep) print }
            ' circ_only.fastq > blat.fastq

            BLAT_READS=\$(awk 'END{print NR/4}' blat.fastq || echo 0)
            echo "BLAT reads: \$BLAT_READS"

            if [ "\$BLAT_READS" -gt 0 ]; then
                echo "Converting BLAT reads to FASTA..."
                awk 'NR%4==1{h=\$0; sub(/^@/,"",h); sub(/ .*/,"",h); print ">"h}
                     NR%4==2{print}' blat.fastq > blat.fa

                echo "Running BLAT..."
                blat circRNAs_blat.fa blat.fa blat.psl ${blat_opts}

                echo "Parsing BLAT PSL to counts..."
                python3 << 'EOFPYTHON'
import os

blat_counts = {}

if os.path.exists("blat.psl"):
    with open("blat.psl") as f:
        for line in f:
            if not line.strip() or line.startswith("psLayout") or line.startswith("match"):
                continue
            cols = line.rstrip().split("\\t")
            if len(cols) < 14:
                continue
            # PSL: tName in col 14 (0-based index 13)
            tName = cols[13]
            blat_counts[tName] = blat_counts.get(tName, 0) + 1

with open("blat_counts.tsv", "w") as out:
    for cid, cnt in blat_counts.items():
        out.write(f"{cid}\\t{cnt}\\n")
EOFPYTHON
            else
                echo "No BLAT reads; leaving blat_counts.tsv empty"
                > blat_counts.tsv
            fi
        else
            echo "No BLAT candidates"
            > blat_counts.tsv
        fi
    else
        echo "BLAT rescue disabled"
        > blat_counts.tsv
    fi

    ########################################
    # 10) Merge BAMs & build debug + ground truth BED12
    ########################################

    if [ -s rescue.bam ]; then
        samtools merge -f merged.bam initial.bam rescue.bam
    else
        cp initial.bam merged.bam
    fi

    python3 << 'EOFPYTHON'
import os

min_reads        = int(${min_reads})
mapq_threshold   = int(${mapq_threshold})
initial_tile     = int(${initial_tile})
rescue_tile      = int(${rescue_tile})
rescue_enabled   = int(${rescue_enabled})
rescue_opts      = "${rescue_opts}"
blat_enabled     = int(${blat_enabled})

orig_bed = os.environ.get("ORIG_BED")
if not orig_bed or not os.path.exists(orig_bed):
    raise SystemExit("ORIG_BED missing")

initial = {}
if os.path.exists("initial_counts.tsv"):
    with open("initial_counts.tsv") as f:
        for line in f:
            if not line.strip():
                continue
            cid, cnt = line.rstrip().split()
            initial[cid] = int(cnt)

rescue = {}
if os.path.exists("rescue_counts.tsv"):
    with open("rescue_counts.tsv") as f:
        for line in f:
            if not line.strip():
                continue
            cid, cnt = line.rstrip().split()
            rescue[cid] = int(cnt)

blat = {}
if os.path.exists("blat_counts.tsv"):
    with open("blat_counts.tsv") as f:
        for line in f:
            if not line.strip():
                continue
            cid, cnt = line.rstrip().split()
            blat[cid] = int(cnt)

headers = {}
if os.path.exists("header_counts.tsv"):
    with open("header_counts.tsv") as f:
        for line in f:
            if not line.strip():
                continue
            key, cnt = line.rstrip().split()
            headers[key] = int(cnt)

expected = {}
if os.path.exists("expected_reads.tsv"):
    with open("expected_reads.tsv") as f:
        next(f, None)
        for line in f:
            if not line.strip():
                continue
            cid, e = line.rstrip().split()
            try:
                expected[cid] = float(e)
            except ValueError:
                expected[cid] = 0.0

def id_to_key(cid: str) -> str:
    k = cid.split("|", 1)[0]
    k = k.split(".", 1)[0]
    return k

total_circs = 0
kept_circs = 0
filtered_circs = 0
rescued_circs = 0
blat_rescued_circs = 0

with open(orig_bed) as bed_in, \\
     open("ground_truth_circular.bed12", "w") as bed_out, \\
     open("ground_truth_counts.debug.tsv", "w") as dbg:

    print("\\t".join([
        "id", "key",
        "init_cnt", "rescue_cnt", "blat_cnt", "total_cnt",
        "fastq_headers", "expected_reads"
    ]), file=dbg)

    for line in bed_in:
        if line.startswith("#") or line.startswith("track"):
            continue
        fields = line.rstrip().split()
        if len(fields) < 4:
            continue

        total_circs += 1
        cid = fields[3]
        key = id_to_key(cid)

        init_cnt = initial.get(cid, 0)
        resc_cnt = rescue.get(cid, 0)
        blat_cnt = blat.get(cid, 0)
        total_cnt = init_cnt + resc_cnt + blat_cnt
        hdr = headers.get(key, 0)
        exp = expected.get(cid, 0.0)

        print("\\t".join([
            cid, key,
            str(init_cnt),
            str(resc_cnt),
            str(blat_cnt),
            str(total_cnt),
            str(hdr),
            f"{exp:.3f}"
        ]), file=dbg)

        if total_cnt >= min_reads:
            # store *raw* total read count as score
            score = total_cnt
            if len(fields) >= 5:
                fields[4] = str(score)
            print("\t".join(fields), file=bed_out)
            kept_circs += 1
            if init_cnt == 0 and resc_cnt > 0:
                rescued_circs += 1
            if init_cnt == 0 and resc_cnt == 0 and blat_cnt > 0:
                blat_rescued_circs += 1
        else:
            filtered_circs += 1

retention = (kept_circs / total_circs * 100.0) if total_circs > 0 else 0.0

with open("filtering_summary.txt", "w") as out:
    out.write(f"Ground Truth Filtering Summary - ${run_name}\\n")
    out.write("="*50 + "\\n")
    out.write("Parameters:\\n")
    out.write(f"  MAPQ threshold (initial): {mapq_threshold}\\n")
    out.write(f"  Initial tiling: {initial_tile}x\\n")
    out.write(f"  Rescue tiling: {rescue_tile}x\\n")
    out.write(f"  Rescue enabled: {rescue_enabled}\\n")
    out.write(f"  Rescue options: {rescue_opts}\\n")
    out.write(f"  BLAT enabled: {blat_enabled}\\n")
    out.write(f"  BLAT tiling: ${blat_tile}x\\n")
    out.write(f"  BLAT opts: ${blat_opts}\\n")
    out.write(f"  Min reads: {min_reads}\\n")
    out.write("="*50 + "\\n")
    out.write(f"Total circRNAs in original BED12: {total_circs}\\n")
    out.write(f"Kept (>= {min_reads} reads): {kept_circs}\\n")
    out.write(f"Filtered (< {min_reads} reads): {filtered_circs}\\n")
    out.write(f"Successfully rescued (minimap2): {rescued_circs}\\n")
    out.write(f"Successfully rescued (BLAT-only): {blat_rescued_circs}\\n")
    out.write(f"Retention rate: {retention:.1f}%\\n")

print(f"Created ground truth: {kept_circs}/{total_circs} ({retention:.1f}%)")
print(f"  Initial nonzero: {len(initial)} circRNAs")
print(f"  Rescue nonzero:  {len(rescue)} circRNAs")
print(f"  BLAT nonzero:    {len(blat)} circRNAs")
print(f"  Rescued (minimap2): {rescued_circs} circRNAs")
print(f"  Rescued (BLAT-only): {blat_rescued_circs} circRNAs")
EOFPYTHON

    ########################################
    # 11) Final BAM for kept circRNAs
    ########################################

    cut -f4 ground_truth_circular.bed12 > kept_ids.txt || > kept_ids.txt

    samtools view -h merged.bam | \\
        awk 'BEGIN{while((getline<"kept_ids.txt")>0) ids[\$0]=1}
             /^@/ || (\$3 in ids)' | \\
        samtools view -bS - | samtools sort -@ ${task.cpus} -o ground_truth.sorted.bam -
    samtools index ground_truth.sorted.bam

    echo "Ground truth BED12 and BAM generated."
    cat filtering_summary.txt
    """
}
