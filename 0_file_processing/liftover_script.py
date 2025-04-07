#!/usr/bin/env python3

import argparse
import os
import re
import subprocess
import sys
from Bio import SeqIO

def check_dependencies():
    """Check if required dependencies are installed."""
    # Check for liftOver
    try:
        subprocess.run(['which', 'liftOver'], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except subprocess.CalledProcessError:
        print("ERROR: liftOver tool not found in PATH.")
        print("Please download it from: https://genome.ucsc.edu/goldenPath/help/hgDownload.html")
        print("And make sure it's in your PATH.")
        return False
    
    return True

def download_chain_file(output_dir):
    """Download the mm9 to mm10 chain file if it doesn't exist."""
    chain_file = os.path.join(output_dir, "mm9ToMm10.over.chain")
    
    if os.path.exists(chain_file):
        print(f"Chain file already exists: {chain_file}")
        return chain_file
    
    # Download gzipped chain file
    gzipped_chain = chain_file + ".gz"
    print(f"Downloading mm9 to mm10 chain file...")
    
    try:
        subprocess.run([
            "wget", 
            "https://hgdownload.soe.ucsc.edu/goldenPath/mm9/liftOver/mm9ToMm10.over.chain.gz",
            "-O", gzipped_chain
        ], check=True)
        
        # Decompress the file
        subprocess.run(["gunzip", "-f", gzipped_chain], check=True)
        print(f"Chain file downloaded and unzipped: {chain_file}")
        return chain_file
    
    except subprocess.CalledProcessError as e:
        print(f"ERROR: Failed to download or unzip chain file: {e}")
        return None

def extract_coordinates(fasta_file, output_dir):
    """Extract mm9 coordinates from circRNA FASTA headers and create BED file."""
    mm9_bed_file = os.path.join(output_dir, "mm9_circrnas.bed")
    
    print(f"Extracting coordinates from FASTA headers: {fasta_file}")
    
    patterns = [
        # Standard format: mmu_circ_0000001|chr1:5114357-5133932+|NM_133826|Atp6v1h
        r'\|([^|:]+):(\d+)-(\d+)([+-])',
        # Alternative format: chr1_5114357-5133932_+
        r'([^_]+)_(\d+)-(\d+)_([+-])',
        # Last resort: just look for chr and coordinates
        r'(chr\w+)[^0-9]*(\d+)[^0-9]*(\d+).*?([+-])'
    ]
    
    coord_count = 0
    with open(mm9_bed_file, "w") as out:
        for record in SeqIO.parse(fasta_file, "fasta"):
            matched = False
            
            for pattern in patterns:
                m = re.search(pattern, record.id)
                if m:
                    chrom, start, end, strand = m.groups()
                    # BED format is 0-based, but circRNA coordinates might be 1-based
                    # Subtract 1 from start to convert to 0-based
                    start_0based = int(start) - 1
                    out.write(f"{chrom}\t{start_0based}\t{end}\t{record.id}\t0\t{strand}\n")
                    coord_count += 1
                    matched = True
                    break
            
            if not matched:
                print(f"WARNING: Could not extract coordinates from header: {record.id}")
    
    print(f"Extracted coordinates for {coord_count} circRNAs")
    return mm9_bed_file if coord_count > 0 else None

def run_liftover(mm9_bed_file, chain_file, output_dir):
    """Run liftOver to convert coordinates from mm9 to mm10."""
    mm10_bed_file = os.path.join(output_dir, "mm10_circrnas.bed")
    unmapped_file = os.path.join(output_dir, "unmapped_circrnas.bed")
    
    print(f"Running liftOver to convert coordinates from mm9 to mm10...")
    
    try:
        subprocess.run([
            "liftOver", 
            mm9_bed_file, 
            chain_file, 
            mm10_bed_file, 
            unmapped_file
        ], check=True)
        
        # Check results
        mapped_count = 0
        with open(mm10_bed_file) as f:
            mapped_count = sum(1 for _ in f)
        
        unmapped_count = 0
        with open(unmapped_file) as f:
            unmapped_count = sum(1 for _ in f if not _.startswith("#"))
        
        print(f"LiftOver results: {mapped_count} mapped, {unmapped_count} unmapped")
        return mm10_bed_file
    
    except subprocess.CalledProcessError as e:
        print(f"ERROR: liftOver failed: {e}")
        return None

def update_fasta(original_fasta, mm10_bed_file, output_dir):
    """Create new FASTA file with updated mm10 coordinates."""
    updated_fasta = os.path.join(output_dir, "mm10_circrnas.fasta")
    
    print(f"Creating updated FASTA file with mm10 coordinates...")
    
    # Load the liftOver results
    mm10_coords = {}
    with open(mm10_bed_file) as f:
        for line in f:
            fields = line.strip().split("\t")
            original_id = fields[3]
            new_chrom = fields[0]
            # Convert back to 1-based for display
            new_start = int(fields[1]) + 1
            new_end = fields[2]
            strand = fields[5]
            mm10_coords[original_id] = (new_chrom, new_start, new_end, strand)
    
    # Update FASTA headers
    records = []
    mapped_count = 0
    for record in SeqIO.parse(original_fasta, "fasta"):
        if record.id in mm10_coords:
            chrom, start, end, strand = mm10_coords[record.id]
            # Update the description to include mm10 coordinates
            record.description = f"{record.id} [mm10:{chrom}:{start}-{end}{strand}]"
            mapped_count += 1
        else:
            record.description = f"{record.id} [UNMAPPED_IN_MM10]"
        records.append(record)
    
    # Write updated FASTA
    with open(updated_fasta, "w") as out:
        SeqIO.write(records, out, "fasta")
    
    print(f"Created updated FASTA with mm10 coordinates: {updated_fasta}")
    print(f"Updated {mapped_count} sequences with mm10 coordinates")
    
    return updated_fasta

def main():
    parser = argparse.ArgumentParser(description='Convert circRNA coordinates from mm9 to mm10')
    parser.add_argument('--fasta', required=True, help='Input circRNA FASTA file with mm9 coordinates')
    parser.add_argument('--output', default='liftover_output', help='Output directory')
    args = parser.parse_args()
    
    # Check if liftOver is installed
    if not check_dependencies():
        sys.exit(1)
    
    # Create output directory if needed
    os.makedirs(args.output, exist_ok=True)
    
    # Download chain file
    chain_file = download_chain_file(args.output)
    if not chain_file:
        sys.exit(1)
    
    # Extract coordinates from FASTA headers
    mm9_bed_file = extract_coordinates(args.fasta, args.output)
    if not mm9_bed_file:
        print("ERROR: Failed to extract coordinates from FASTA headers")
        sys.exit(1)
    
    # Run liftOver
    mm10_bed_file = run_liftover(mm9_bed_file, chain_file, args.output)
    if not mm10_bed_file:
        sys.exit(1)
    
    # Create updated FASTA
    updated_fasta = update_fasta(args.fasta, mm10_bed_file, args.output)
    
    print("Conversion complete!")
    print(f"Updated FASTA file: {updated_fasta}")
    print(f"You can now use this file with your mm10 GTF annotation.")
    
    # Return the path to the updated FASTA for use in other scripts
    return updated_fasta

if __name__ == "__main__":
    main()