#!/usr/bin/env python3
"""
Improved circRNA Simulator with Biological Features

This script simulates circRNAs with biologically accurate features including:
- Proper distribution of circRNA types (eciRNA, EIciRNA, ciRNA, intergenic)
- Accurate mature length distributions by type
- Biologically realistic splice site patterns by extracting actual splice sites from the genome
- Rolling circle parameters (copy numbers) based on biological data

Features are based on real data analysis from database extraction.
"""

import sys
import os
import numpy as np
import pandas as pd
from tqdm import tqdm
from collections import defaultdict
import random
import logging
import argparse
import subprocess
import tempfile
import shutil

# Function to set up logging
def setup_logging(output_dir):
    os.makedirs(output_dir, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(os.path.join(output_dir, 'simulation.log')),
            logging.StreamHandler()
        ]
    )

# Function to reverse complement a DNA sequence
def reverse_complement(seq):
    """Return reverse complement of a sequence."""
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 
                  'a': 't', 'c': 'g', 'g': 'c', 't': 'a',
                  'N': 'N', 'n': 'n'}
    return ''.join(complement.get(base, base) for base in reversed(seq))

# Function to parse FASTA files using a generator
def parse_fasta(fasta_file):
    """Parse a FASTA file and yield sequence records."""
    with open(fasta_file, 'r') as f:
        header = None
        sequence = []
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if header:
                    yield header, ''.join(sequence)
                header = line[1:]
                sequence = []
            else:
                sequence.append(line)
        if header:
            yield header, ''.join(sequence)

# Class to load and access genome sequences
class GenomeLoader:
    def __init__(self, genome_file):
        """Load genome sequences from FASTA file into memory."""
        self.sequences = {}
        logging.info(f"Loading genome from {genome_file}...")
        for header, sequence in parse_fasta(genome_file):
            # Extract chromosome name from header
            chrom = header.split()[0]
            self.sequences[chrom] = sequence
        logging.info(f"Loaded {len(self.sequences)} chromosomes/scaffolds")
    
    def get_sequence(self, chrom, start, end, strand='+'):
        """Get sequence from genome with 0-based coordinates (start, end)."""
        if chrom not in self.sequences:
            logging.warning(f"Chromosome {chrom} not found in genome")
            return ""
        
        # Adjust coordinates (make sure they're within bounds)
        start = max(0, start)
        end = min(len(self.sequences[chrom]), end)
        
        # Get sequence
        seq = self.sequences[chrom][start:end]
        
        # Reverse complement if needed
        if strand == '-':
            seq = reverse_complement(seq)
        
        return seq

# Class to parse and organize GTF data
class GTFParser:
    def __init__(self, gtf_file, bedtools_path=None):
        """Initialize with GTF file path."""
        self.gtf_file = gtf_file
        self.bedtools_path = bedtools_path or 'bedtools'
        self.temp_dir = None
        self.genes_gtf = None
        self.exons_gtf = None
        
        # Parse the GTF file and prepare data structures
        self.prepare_gtf_data()
    
    def prepare_gtf_data(self):
        """Extract and organize gene, exon, and intron features."""
        # Create temporary directory
        self.temp_dir = tempfile.mkdtemp()
        logging.info(f"Preparing GTF data in temporary directory: {self.temp_dir}")
        
        # Extract genes and exons
        self.genes_gtf = os.path.join(self.temp_dir, "genes.gtf")
        self.exons_gtf = os.path.join(self.temp_dir, "exons.gtf")
        
        # Use grep to extract genes and exons (faster than parsing the whole GTF)
        with open(self.gtf_file, 'r') as gtf, \
             open(self.genes_gtf, 'w') as genes, \
             open(self.exons_gtf, 'w') as exons:
            
            for line in gtf:
                if line.startswith('#'):
                    continue
                if '\tgene\t' in line:
                    genes.write(line)
                elif '\texon\t' in line:
                    exons.write(line)
        
        logging.info("Extracted gene and exon features from GTF")
        
        # Load transcript structure
        self.tx_structure = self.load_transcript_structure()
        logging.info(f"Loaded transcript structure for {len(self.tx_structure)} genes")
    
    def load_transcript_structure(self):
        """Load transcript structure from exon GTF."""
        tx_structure = defaultdict(dict)
        
        with open(self.exons_gtf, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                
                fields = line.strip().split('\t')
                if len(fields) < 9:
                    continue
                
                # Extract basic information
                chrom = fields[0]
                feature_type = fields[2]
                start = int(fields[3])
                end = int(fields[4])
                strand = fields[6]
                attributes = fields[8]
                
                # Skip non-exon features
                if feature_type != 'exon':
                    continue
                
                # Parse attributes
                attr_dict = {}
                for attr in attributes.split(';'):
                    attr = attr.strip()
                    if not attr:
                        continue
                    
                    if '"' in attr:
                        key, value = attr.split('"', 1)
                        key = key.strip()
                        value = value.strip('"').strip()
                    else:
                        parts = attr.split()
                        if len(parts) < 2:
                            continue
                        key, value = parts[0], parts[1]
                        value = value.strip('"')
                    
                    attr_dict[key] = value
                
                # Extract IDs
                gene_id = attr_dict.get('gene_id', '')
                transcript_id = attr_dict.get('transcript_id', '')
                gene_type = attr_dict.get('gene_type', attr_dict.get('gene_biotype', ''))
                
                # Skip if missing essential IDs
                if not gene_id or not transcript_id:
                    continue
                
                # Create exon entry
                exon_entry = {
                    'chrom': chrom,
                    'start': start,
                    'end': end,
                    'strand': strand,
                    'gene_id': gene_id,
                    'transcript_id': transcript_id,
                    'gene_type': gene_type
                }
                
                # Store in transcript structure
                if transcript_id not in tx_structure[gene_id]:
                    tx_structure[gene_id][transcript_id] = []
                tx_structure[gene_id][transcript_id].append(exon_entry)
        
        # Sort exons by position for each transcript
        for gene_id in tx_structure:
            for tx_id in tx_structure[gene_id]:
                tx_structure[gene_id][tx_id].sort(key=lambda x: x['start'])
        
        return tx_structure
    
    def is_protein_coding(self, gene_id):
        """Check if a gene is protein-coding."""
        if gene_id in self.tx_structure:
            for tx_id in self.tx_structure[gene_id]:
                if self.tx_structure[gene_id][tx_id]:
                    gene_type = self.tx_structure[gene_id][tx_id][0].get('gene_type', '')
                    return gene_type == 'protein_coding'
        return False
    
    def get_genes_with_sufficient_exons(self, min_exons=2, protein_coding_only=True):
        """Get list of genes with at least the specified number of exons."""
        valid_genes = []
        
        for gene_id in self.tx_structure:
            # Skip non-protein-coding if specified
            if protein_coding_only and not self.is_protein_coding(gene_id):
                continue
            
            # Check if any transcript has sufficient exons
            for tx_id, exons in self.tx_structure[gene_id].items():
                if len(exons) >= min_exons:
                    valid_genes.append((gene_id, tx_id))
                    break
        
        return valid_genes
    
    def classify_region(self, chrom, start, end, strand):
        """
        Classify a genomic region as exonic, intronic, or intergenic.
        
        Using bedtools logic:
        - If overlaps fully with exons (-s -split -f 1.0) -> exonic (eciRNA)
        - If overlaps partially with exons (-s -split) -> exon-intron (EIciRNA)
        - If overlaps with genes but not exons (-s) -> intronic (ciRNA)
        - If doesn't overlap with genes -> intergenic
        """
        # Create a temporary BED file for the region
        region_bed = os.path.join(self.temp_dir, "temp_region.bed")
        with open(region_bed, 'w') as f:
            f.write(f"{chrom}\t{start}\t{end}\t.\t.\t{strand}\n")
        
        # Check gene overlap (genic vs intergenic)
        gene_overlap_file = os.path.join(self.temp_dir, "gene_overlap.bed")
        gene_cmd = f"{self.bedtools_path} intersect -a {region_bed} -b {self.genes_gtf} -s -wa -wb > {gene_overlap_file}"
        subprocess.run(gene_cmd, shell=True, check=False)
        
        # Check if region overlaps with any genes
        with open(gene_overlap_file, 'r') as f:
            gene_overlaps = f.readlines()
        
        if not gene_overlaps:
            return "intergenic"
        
        # Check exon overlap with -split (allows for partial overlap)
        exon_split_file = os.path.join(self.temp_dir, "exon_split_overlap.bed")
        exon_split_cmd = f"{self.bedtools_path} intersect -a {region_bed} -b {self.exons_gtf} -s -split -wa -wb > {exon_split_file}"
        subprocess.run(exon_split_cmd, shell=True, check=False)
        
        # Check if region overlaps with any exons
        with open(exon_split_file, 'r') as f:
            exon_split_overlaps = f.readlines()
        
        if not exon_split_overlaps:
            return "ciRNA"  # Intronic circRNA
        
        # Check for complete exon overlap
        full_exon_file = os.path.join(self.temp_dir, "full_exon_overlap.bed")
        full_exon_cmd = f"{self.bedtools_path} intersect -a {region_bed} -b {self.exons_gtf} -s -split -f 1.0 -wa > {full_exon_file}"
        subprocess.run(full_exon_cmd, shell=True, check=False)
        
        # Check if region fully overlaps with exons
        with open(full_exon_file, 'r') as f:
            full_exon_overlaps = f.readlines()
        
        if full_exon_overlaps:
            return "eciRNA"  # Exonic circRNA
        else:
            return "EIciRNA"  # Exon-intron circRNA
    
    def extract_splice_sites(self, genome, chrom, start, end, strand):
        """
        Extract splice site dinucleotides from genome reference.
        For circRNAs, we need to find the donor and acceptor sites at the backsplice junction.
        Returns: (donor_site, acceptor_site, splice_type)
        """
        # For circRNAs, the donor site is at the end coordinate and the acceptor is at the start
        if strand == '+':
            # Get 2bp upstream and 2bp downstream of junction
            donor_region = (end - 2, end + 2)  # GT at the end
            acceptor_region = (start - 2, start + 2)  # AG at the start
        else:
            # For negative strand, swap donor and acceptor roles but keep same coordinates
            donor_region = (start - 2, start + 2)  # AC at the start (reverse complement of GT)
            acceptor_region = (end - 2, end + 2)  # CT at the end (reverse complement of AG)
        
        # Get the sequences
        donor_seq = genome.get_sequence(chrom, donor_region[0], donor_region[1], '+')
        acceptor_seq = genome.get_sequence(chrom, acceptor_region[0], acceptor_region[1], '+')
        
        # Extract the dinucleotide splice sites
        if strand == '+':
            donor_dinuc = donor_seq[2:4]  # Get the GT
            acceptor_dinuc = acceptor_seq[0:2]  # Get the AG
        else:
            # For negative strand, get the reverse complement
            donor_dinuc = reverse_complement(acceptor_seq[0:2])  # Get GT (RC of AC)
            acceptor_dinuc = reverse_complement(donor_seq[2:4])  # Get AG (RC of CT)
        
        # Classify splice site type
        if (donor_dinuc, acceptor_dinuc) == ('GT', 'AG'):
            splice_type = "GT-AG"
        elif (donor_dinuc, acceptor_dinuc) == ('GC', 'AG'):
            splice_type = "GC-AG"
        elif (donor_dinuc, acceptor_dinuc) == ('AT', 'AC'):
            splice_type = "AT-AC"
        else:
            splice_type = "Non-canonical"
        
        return donor_dinuc, acceptor_dinuc, splice_type
    
    def create_suitable_splice_sites(self, genome, chrom, start, end, strand, desired_type="GT-AG"):
        """
        Find positions near the specified start/end that would create the desired splice site type.
        This is used to adjust the circRNA boundaries to ensure canonical splice sites.
        Returns: (new_start, new_end)
        """
        # Define search range (search up to 10bp in each direction)
        max_adjustment = 10
        best_positions = None
        
        # Try positions within the search range
        for start_adj in range(-max_adjustment, max_adjustment + 1):
            for end_adj in range(-max_adjustment, max_adjustment + 1):
                new_start = start + start_adj
                new_end = end + end_adj
                
                # Get splice sites at these positions
                donor, acceptor, splice_type = self.extract_splice_sites(genome, chrom, new_start, new_end, strand)
                
                # If we found the desired type, return the adjusted positions
                if splice_type == desired_type:
                    return new_start, new_end
        
        # If no exact match found, return original positions
        return start, end
    
    def cleanup(self):
        """Clean up temporary files."""
        if self.temp_dir and os.path.exists(self.temp_dir):
            shutil.rmtree(self.temp_dir)
            self.temp_dir = None

# Class to simulate circRNAs
class CircRNASimulator:
    def __init__(self, gtf_parser, genome_loader, output_dir, bedtools_path=None):
        """Initialize simulator with parsed GTF data and genome."""
        self.gtf_parser = gtf_parser
        self.genome = genome_loader
        self.output_dir = output_dir
        self.bedtools_path = bedtools_path or 'bedtools'
        
        # Set up biological feature distributions
        self.setup_biological_features()
    
    def setup_biological_features(self):
        """Set up realistic biological feature distributions based on analyzed data."""
        # circRNA type distribution (based on the visualization_report.html data)
        self.type_distribution = {
            'eciRNA': 0.545,    # Exonic circRNAs (54.5%)
            'EIciRNA': 0.370,   # Exon-intron circRNAs (37.0%)
            'ciRNA': 0.014,     # Intronic circRNAs (1.4%)
            'intergenic': 0.071 # Intergenic circRNAs (7.1%)
        }
        
        # Length distributions by type (based on visualization_report.html data)
        self.length_distributions = {
            'eciRNA': {'median': 393, 'q1': 200, 'q3': 700, 'min': 80, 'max': 3000},
            'EIciRNA': {'median': 519.5, 'q1': 300, 'q3': 900, 'min': 100, 'max': 4000},
            'ciRNA': {'median': 314.5, 'q1': 150, 'q3': 600, 'min': 80, 'max': 2000},
            'intergenic': {'median': 365, 'q1': 180, 'q3': 650, 'min': 80, 'max': 2500}
        }
        
        # Splice site distributions (based on real data)
        self.splice_site_distribution = {
            'GT-AG': 0.90,  # Canonical (majority)
            'GC-AG': 0.05,  # Non-canonical but common
            'AT-AC': 0.01,  # U12-type
            'non-canonical': 0.04  # Truly non-canonical
        }
        
        # Update exon count distribution with precise percentages
        self.exon_count_distribution = {
            'eciRNA': {
                1: 0.115,   # 11.5% single exon
                2: 0.323,   # 32.3% two exons
                3: 0.445,   # 44.5% three exons
                4: 0.149,   # 14.9% four exons
                '5+': 0.193 # 19.3% five or more exons
            },
            'EIciRNA': {
                1: 0.082,   # 8.2% single exon
                2: 0.179,   # 17.9% two exons
                3: 0.213,   # 21.3% three exons
                4: 0.183,   # 18.3% four exons
                '5+': 0.343 # 34.3% five or more exons
            },
            'ciRNA': {
                1: 0.594,   # 59.4% single exon
                2: 0.241,   # 24.1% two exons
                3: 0.094,   # 9.4% three exons
                4: 0.0,     # 0% four exons
                '5+': 0.0   # 0% five or more exons
            },
            'intergenic': {
                1: 0.350,   # 35.0% single exon
                2: 0.289,   # 28.9% two exons
                3: 0.197,   # 19.7% three exons
                4: 0.093,   # 9.3% four exons
                '5+': 0.069 # 6.9% five or more exons
            }
        }
        
        # Rolling circle parameters (from simulation_parameters.txt)
        # Create a proper distribution of copy numbers
        self.mean_copy_number = 13.50  # from simulation_parameters.txt
        
        # Create a copy number distribution from the data in simulation_parameters.txt
        self.copy_number_distribution = self.generate_copy_number_distribution()
    
    def generate_copy_number_distribution(self):
        """Generate a distribution of copy numbers based on the real data."""
        # Key copy number ranges
        ranges = [
            (1, 3),    # Low copy numbers (1-3)
            (3, 6),    # Medium-low copy numbers (3-6)
            (6, 10),   # Medium copy numbers (6-10)
            (10, 15),  # Medium-high copy numbers (10-15)
            (15, 20),  # High copy numbers (15-20)
            (20, 50),  # Very high copy numbers (20-50)
            (50, 100), # Extreme copy numbers (50-100)
            (100, 1500) # Rare extremely high copy numbers (>100)
        ]
        
        # Probabilities based on the real data distribution
        probabilities = [0.30, 0.30, 0.20, 0.10, 0.05, 0.03, 0.015, 0.005]
        
        return {"ranges": ranges, "probabilities": probabilities}
    
    def sample_from_distribution(self, distribution):
        """
        Sample a value from a discrete probability distribution.
        
        Ensures probabilities sum to 1 by normalizing.
        
        Args:
            distribution (dict): Dictionary of items and their probabilities
        
        Returns:
            The selected item
        """
        items = list(distribution.keys())
        probabilities = list(distribution.values())
        
        # Normalize probabilities to ensure they sum to 1
        total = sum(probabilities)
        normalized_probabilities = [p / total for p in probabilities]
        
        return np.random.choice(items, p=normalized_probabilities)
    
    def generate_copy_number(self):
        """Generate a realistic copy number based on the distribution."""
        # Select a range
        idx = np.random.choice(len(self.copy_number_distribution["ranges"]), 
                              p=self.copy_number_distribution["probabilities"])
        
        # Get the selected range
        min_val, max_val = self.copy_number_distribution["ranges"][idx]
        
        # Generate a value within that range
        if idx == len(self.copy_number_distribution["ranges"]) - 1:
            # For the highest range, use a long-tailed distribution
            # This creates a small chance of extremely high copy numbers
            base = np.random.uniform(min_val, min_val + 40)  # Most values near the beginning
            if np.random.random() < 0.1:  # 10% chance of a very high value
                return int(base + np.random.exponential(scale=100))
            return int(base)
        else:
            # For other ranges, uniform distribution
            return int(np.random.uniform(min_val, max_val))
    
    def generate_block_structure(self, circ):
        if circ['type'] == 'EIciRNA':
            # EIciRNA block distribution
            block_choice = np.random.choice(
                [1, 2, 3, 4, 5], 
                p=[0.082, 0.179, 0.213, 0.183, 0.343]
            )
            
            # Rest of the block generation logic remains the same
            total_length = circ['mature_length']
            
            if block_choice == 1:
                return [total_length], [0]
        
        # Multiple blocks
        block_sizes = [total_length // block_choice + random.randint(-30, 30) for _ in range(block_choice-1)]
        # Adjust last block to ensure total length matches
        block_sizes.append(total_length - sum(block_sizes))
        
        # Calculate block starts
        block_starts = [sum(block_sizes[:i]) for i in range(block_choice)]
        
        return block_sizes, block_starts

    def generate_length(self, circ_type):
        """Generate a random length for a given circRNA type following the real distribution."""
        params = self.length_distributions[circ_type]
        
        # Generate from a skewed normal-like distribution based on quartiles
        # This creates a more realistic distribution than pure normal
        r = np.random.random()
        if r < 0.25:
            # Q1 range
            return int(np.random.uniform(params['min'], params['q1']))
        elif r < 0.5:
            # Q1-median range
            return int(np.random.uniform(params['q1'], params['median']))
        elif r < 0.75:
            # median-Q3 range
            return int(np.random.uniform(params['median'], params['q3']))
        else:
            # Above Q3 (long tail)
            if np.random.random() < 0.9:
                # 90% in normal Q3-max range
                return int(np.random.uniform(params['q3'], params['q3']*2))
            else:
                # 10% in extended tail
                return int(np.random.uniform(params['q3']*2, params['max']))
    
    def generate_exon_count(self, circ_type=None):
        """
        Generate number of exons for circRNAs.
        
        Args:
            circ_type (str, optional): Type of circRNA 
            (eciRNA, EIciRNA, ciRNA, intergenic)
        
        Returns:
            int: Number of exons
        """
        # If no type specified, use a default or raise an error
        if circ_type is None:
            raise ValueError("CircRNA type must be specified")
        
        # Get distribution for the specific circRNA type
        dist = self.exon_count_distribution.get(circ_type, {})
        
        # If no distribution found, raise an error
        if not dist:
            raise ValueError(f"No exon count distribution found for type: {circ_type}")
        
        # Sample from the distribution
        count = self.sample_from_distribution(dist)
        
        # If '5+', generate a random number between 5 and 10
        if count == '5+':
            return np.random.randint(5, 11)
        
        return int(count)  # Ensure we return an integer
    
    def select_splice_site_by_type(self, circ_type):
        """
        Select splice site type based on the distribution from the visualization
        
        Args:
            circ_type (str): Type of circRNA (eciRNA, EIciRNA, ciRNA, intergenic)
        
        Returns:
            str: Selected splice site type
        """
        # Distributions based on the provided visualization
        type_distributions = {
            'eciRNA': {'GT-AG': 1.0, 'non-canonical': 0.0},
            'EIciRNA': {'GT-AG': 0.998, 'non-canonical': 0.002},
            'ciRNA': {'GT-AG': 0.912, 'non-canonical': 0.088},
            'intergenic': {'GT-AG': 0.505, 'non-canonical': 0.495}
        }
        
        # Get distribution for the specific circRNA type
        dist = type_distributions.get(circ_type, {'GT-AG': 1.0, 'non-canonical': 0.0})
        
        # Select splice site type based on distribution
        r = np.random.random()
        cumulative = 0
        for splice_type, prob in dist.items():
            cumulative += prob
            if r < cumulative:
                return splice_type
        
        # Fallback to GT-AG if something goes wrong
        return 'GT-AG'
    
    def select_non_canonical_type(self):
        """
        Selects a specific non-canonical splice site type.
        
        Returns:
            str: A non-canonical splice site type
        """
        # Distribution of non-canonical types
        non_canonical_types = {
            'GC-AG': 0.5,     # Most common non-canonical
            'AT-AC': 0.1,     # Less common U12-type
            'other': 0.4       # Truly rare non-canonical
        }
        
        r = np.random.random()
        cumulative = 0
        for splice_type, prob in non_canonical_types.items():
            cumulative += prob
            if r < cumulative:
                return splice_type
        
        return 'GC-AG'  # Default fallback

    def simulate_eciRNA(self, num_circrnas, genome, gtf_parser):
        """Simulate a batch of exonic circRNAs (entirely from exons)."""
        logging.info(f"Simulating {num_circrnas} eciRNA (exonic) circRNAs...")
        
        # Get genes with sufficient exons
        valid_genes = gtf_parser.get_genes_with_sufficient_exons(min_exons=2, protein_coding_only=True)
        if not valid_genes:
            logging.error("No valid genes found with sufficient exons")
            return []
        
        results = []
        for _ in tqdm(range(num_circrnas)):
            # Generate exon count for eciRNA
            exon_count = self.generate_exon_count('eciRNA')
            
            # Randomly select a gene and transcript
            gene_id, tx_id = random.choice(valid_genes)
            exons = gtf_parser.tx_structure[gene_id][tx_id]
            
            # Ensure we have enough exons for the desired count
            if len(exons) < exon_count:
                continue
            
            # Randomly select start exon index to accommodate exon count
            max_start_idx = len(exons) - exon_count
            start_idx = np.random.randint(0, max_start_idx + 1)
            end_idx = start_idx + exon_count - 1
            
            # Get the selected exons
            selected_exons = exons[start_idx:end_idx+1]
            
            # Get basic information
            chrom = selected_exons[0]['chrom']
            strand = selected_exons[0]['strand']
            circ_start = selected_exons[0]['start']
            circ_end = selected_exons[-1]['end']
            
            # Select splice site type based on eciRNA distribution
            desired_splice_type = self.select_splice_site_by_type('eciRNA')
            
            # Try to find suitable splice sites
            circ_start, circ_end = gtf_parser.create_suitable_splice_sites(
                genome, chrom, circ_start, circ_end, strand, desired_splice_type)
            
            # Verify the region is actually eciRNA
            region_type = gtf_parser.classify_region(chrom, circ_start, circ_end, strand)
            if region_type != "eciRNA":
                logging.debug(f"Region classified as {region_type}, not eciRNA, retrying...")
                continue
            
            # Generate the sequence
            circ_seq = ""
            for exon in selected_exons:
                exon_seq = genome.get_sequence(exon['chrom'], exon['start']-1, exon['end'], exon['strand'])
                circ_seq += exon_seq
            
            # Verify if sequence was obtained
            if not circ_seq:
                logging.warning(f"Failed to obtain sequence for eciRNA at {chrom}:{circ_start}-{circ_end}")
                continue
            
            # Generate random copy number for rolling circle structure
            copy_number = self.generate_copy_number()
            
            # Create CircRNA entry - preserve original Ensembl IDs
            # Format: ENSTXXXX|in_ENS|eciRNA|chrom:start-end|strand|length
            circ_id = f"{tx_id}|in_ENS|eciRNA|{chrom}:{circ_start}-{circ_end}|{strand}|{len(circ_seq)}"
            
            # Create pseudo-circular sequence with multiple copies
            pseudo_seq = self.generate_rolling_circle_sequence(circ_seq, copy_number)
            
            # Random expression level (using beta distribution for skewed values)
            tpm = np.random.beta(0.5, 2.0) * 100
            
            results.append({
                'id': circ_id,
                'chrom': chrom,
                'start': circ_start,
                'end': circ_end,
                'strand': strand,
                'type': 'eciRNA',
                'gene_id': gene_id,
                'tx_id': tx_id,
                'circ_seq': circ_seq,
                'pseudo_seq': pseudo_seq,
                'mature_length': len(circ_seq),
                'exon_count': exon_count,
                'copy_number': copy_number,
                'tpm': tpm,
                'splice_type': desired_splice_type
            })
        
        return results
    
    def simulate_EIciRNA(self, num_circrnas, genome, gtf_parser):
        """Simulate exon-intron circRNAs (partially overlapping exons)."""
        logging.info(f"Simulating {num_circrnas} EIciRNA (exon-intron) circRNAs...")
        
        # Get genes with sufficient exons
        valid_genes = gtf_parser.get_genes_with_sufficient_exons(min_exons=2, protein_coding_only=True)
        if not valid_genes:
            logging.error("No valid genes found with sufficient exons")
            return []
        
        results = []
        for _ in tqdm(range(num_circrnas)):
            # Generate exon count for EIciRNA
            exon_count = self.generate_exon_count('EIciRNA')
            
            # Randomly select a gene and transcript
            gene_id, tx_id = random.choice(valid_genes)
            exons = gtf_parser.tx_structure[gene_id][tx_id]
            
            # Ensure we have enough exons for the desired count
            if len(exons) < exon_count:
                continue
            
            # Randomly select start exon index to accommodate exon count
            max_start_idx = len(exons) - exon_count
            start_idx = np.random.randint(0, max_start_idx + 1)
            end_idx = start_idx + exon_count - 1
            
            # Get the selected exons and surrounding introns
            selected_segments = []
            for i in range(start_idx, end_idx + 1):
                # Add exon
                selected_segments.append(exons[i])
                
                # Add intron if not the last exon
                if i < end_idx:
                    next_exon = exons[i + 1]
                    intron_start = exons[i]['end'] + 1
                    intron_end = next_exon['start'] - 1
                    
                    # Optionally skip or partially include intron based on probability
                    if random.random() > 0.3:  # 70% chance of including intron
                        intron_length = intron_end - intron_start + 1
                        intron_seq = genome.get_sequence(
                            exons[i]['chrom'], 
                            intron_start - 1, 
                            intron_end, 
                            exons[i]['strand']
                        )
                        selected_segments.append({
                            'type': 'intron',
                            'start': intron_start,
                            'end': intron_end,
                            'sequence': intron_seq
                        })
            
            # Generate sequence from selected segments
            circ_seq = ""
            for segment in selected_segments:
                if isinstance(segment, dict) and 'sequence' in segment:
                    circ_seq += segment['sequence']
                else:
                    seg_seq = genome.get_sequence(
                        segment['chrom'], 
                        segment['start'] - 1, 
                        segment['end'], 
                        segment['strand']
                    )
                    circ_seq += seg_seq
            
            # Select splice site type based on EIciRNA distribution
            desired_splice_type = self.select_splice_site_by_type('EIciRNA')
            
            # If non-canonical is selected, choose a specific type
            if desired_splice_type == 'non-canonical':
                desired_splice_type = self.select_non_canonical_type()
            
            # Try to find suitable splice sites
            circ_start = selected_segments[0]['start']
            circ_end = selected_segments[-1]['end']
            chrom = selected_segments[0]['chrom']
            strand = selected_segments[0]['strand']
            
            circ_start, circ_end = gtf_parser.create_suitable_splice_sites(
                genome, chrom, circ_start, circ_end, strand, desired_splice_type)
            
            # Verify the region is actually EIciRNA
            region_type = gtf_parser.classify_region(chrom, circ_start, circ_end, strand)
            if region_type != "EIciRNA":
                logging.debug(f"Region classified as {region_type}, not EIciRNA, retrying...")
                continue
            
            # Verify if sequence was obtained
            if not circ_seq:
                logging.warning(f"Failed to obtain sequence for EIciRNA at {chrom}:{circ_start}-{circ_end}")
                continue
            
            # Generate random copy number for rolling circle structure
            copy_number = self.generate_copy_number()
            
            # Create CircRNA entry - preserve original Ensembl IDs
            # Format: ENSTXXXX|in_ENS|EIciRNA|chrom:start-end|strand|length
            circ_id = f"{tx_id}|in_ENS|EIciRNA|{chrom}:{circ_start}-{circ_end}|{strand}|{len(circ_seq)}"
            
            # Create pseudo-circular sequence with multiple copies
            pseudo_seq = self.generate_rolling_circle_sequence(circ_seq, copy_number)
            
            # Random expression level (using beta distribution for skewed values)
            tpm = np.random.beta(0.5, 2.0) * 80  # Slightly lower expression than eciRNA
            
            results.append({
                'id': circ_id,
                'chrom': chrom,
                'start': circ_start,
                'end': circ_end,
                'strand': strand,
                'type': 'EIciRNA',
                'gene_id': gene_id,
                'tx_id': tx_id,
                'circ_seq': circ_seq,
                'pseudo_seq': pseudo_seq,
                'mature_length': len(circ_seq),
                'exon_count': exon_count,
                'copy_number': copy_number,
                'tpm': tpm,
                'splice_type': desired_splice_type
            })
        
        return results
    
    def simulate_ciRNA(self, num_circrnas, genome, gtf_parser):
        """Simulate intronic circRNAs (entirely within introns)."""
        logging.info(f"Simulating {num_circrnas} ciRNA (intronic) circRNAs...")
        
        # Get genes with sufficient exons to have introns
        valid_genes = gtf_parser.get_genes_with_sufficient_exons(min_exons=2, protein_coding_only=True)
        if not valid_genes:
            logging.error("No valid genes found with sufficient exons")
            return []
        
        results = []
        for _ in tqdm(range(num_circrnas)):
            # Generate exon count for ciRNA
            exon_count = self.generate_exon_count('ciRNA')
            
            # Randomly select a gene and transcript
            gene_id, tx_id = random.choice(valid_genes)
            exons = gtf_parser.tx_structure[gene_id][tx_id]
            
            # Ensure we have enough exons to have meaningful introns
            if len(exons) < 2:
                continue
            
            # Select intron(s) based on exon count
            selected_introns = []
            num_introns_to_select = min(exon_count, len(exons) - 1)
            
            # Randomly select intron indices
            intron_indices = random.sample(range(len(exons) - 1), num_introns_to_select)
            
            for idx in intron_indices:
                first_exon = exons[idx]
                second_exon = exons[idx + 1]
                
                # Define the intron region
                if first_exon['end'] < second_exon['start']:
                    intron_start = first_exon['end'] + 1
                    intron_end = second_exon['start'] - 1
                else:
                    intron_start = second_exon['end'] + 1
                    intron_end = first_exon['start'] - 1
                
                # Ensure intron is long enough
                if intron_end - intron_start + 1 < 100:
                    continue
                
                selected_introns.append({
                    'start': intron_start,
                    'end': intron_end,
                    'chrom': first_exon['chrom'],
                    'strand': first_exon['strand']
                })
            
            # If no valid introns found, skip
            if not selected_introns:
                continue
            
            # Choose an intron or combine multiple introns
            if len(selected_introns) == 1:
                intron = selected_introns[0]
                circ_start = intron['start']
                circ_end = intron['end']
            else:
                # Combine multiple introns with some randomness
                circ_start = min(intron['start'] for intron in selected_introns)
                circ_end = max(intron['end'] for intron in selected_introns)
            
            # Generate the sequence
            circ_seq = genome.get_sequence(
                selected_introns[0]['chrom'], 
                circ_start-1, 
                circ_end, 
                selected_introns[0]['strand']
            )
            
            # Verify if sequence was obtained
            if not circ_seq:
                logging.warning(f"Failed to obtain sequence for ciRNA at {selected_introns[0]['chrom']}:{circ_start}-{circ_end}")
                continue
            
            # Generate random copy number for rolling circle structure
            copy_number = self.generate_copy_number()
            
            # Select splice site type based on ciRNA distribution
            desired_splice_type = self.select_splice_site_by_type('ciRNA')
            
            # If non-canonical is selected, choose a specific type
            if desired_splice_type == 'non-canonical':
                desired_splice_type = self.select_non_canonical_type()
            
            # Try to find suitable splice sites
            circ_start, circ_end = gtf_parser.create_suitable_splice_sites(
                genome, selected_introns[0]['chrom'], circ_start, circ_end, selected_introns[0]['strand'], desired_splice_type)
            
            # Verify the region is actually ciRNA
            region_type = gtf_parser.classify_region(selected_introns[0]['chrom'], circ_start, circ_end, selected_introns[0]['strand'])
            if region_type != "ciRNA":
                logging.debug(f"Region classified as {region_type}, not ciRNA, retrying...")
                continue
            
            # Random expression level (using beta distribution for skewed values)
            tpm = np.random.beta(0.3, 2.5) * 60  # Lower expression than eciRNA/EIciRNA
            
            # Create CircRNA entry
            circ_id = f"{tx_id}|in_ENS|ciRNA|{selected_introns[0]['chrom']}:{circ_start}-{circ_end}|{selected_introns[0]['strand']}|{len(circ_seq)}"
            pseudo_seq = self.generate_rolling_circle_sequence(circ_seq, copy_number)
            
            results.append({
                'id': circ_id,
                'chrom': selected_introns[0]['chrom'],
                'start': circ_start,
                'end': circ_end,
                'strand': selected_introns[0]['strand'],
                'type': 'ciRNA',
                'gene_id': gene_id,
                'tx_id': tx_id,
                'circ_seq': circ_seq,
                'pseudo_seq': pseudo_seq,
                'mature_length': len(circ_seq),
                'exon_count': exon_count,
                'copy_number': copy_number,
                'tpm': tpm,
                'splice_type': desired_splice_type
            })
        
        return results
    
    def simulate_intergenic_circRNA(self, num_circrnas, genome, gtf_parser):
        """Simulate intergenic circRNAs (outside of genes)."""
        logging.info(f"Simulating {num_circrnas} intergenic circRNAs...")
        
        # Get chromosomes and their lengths
        chroms = list(genome.sequences.keys())
        
        results = []
        attempts = 0
        max_attempts = num_circrnas * 10  # Limit attempts to avoid infinite loop
        
        while len(results) < num_circrnas and attempts < max_attempts:
            attempts += 1
            
            # Generate exon count for intergenic circRNA
            exon_count = self.generate_exon_count('intergenic')
            
            # Select splice site type based on intergenic distribution
            splice_type_options = {
                'GT-AG': 0.505,
                'non-canonical': 0.495
            }
            splice_type = self.sample_from_distribution(splice_type_options)
            
            # If non-canonical is selected, choose a specific type
            if splice_type == 'non-canonical':
                non_canonical_types = {
                    'GC-AG': 0.5,     # Most common non-canonical
                    'AT-AC': 0.1,     # Less common U12-type
                    'other': 0.4       # Truly rare non-canonical
                }
                splice_type = self.sample_from_distribution(non_canonical_types)
            
            # Randomly select a chromosome
            chrom = random.choice(chroms)
            chrom_length = len(genome.sequences[chrom])
            
            # Skip very small chromosomes
            if chrom_length < 10000:
                continue
            
            # Generate a random length
            circ_length = self.generate_length('intergenic')
            
            # Random start position
            circ_start = np.random.randint(1, chrom_length - circ_length - 1)
            circ_end = circ_start + circ_length - 1
            
            # Random strand
            strand = random.choice(['+', '-'])
            
            # Verify the region is actually intergenic
            region_type = gtf_parser.classify_region(chrom, circ_start, circ_end, strand)
            if region_type != "intergenic":
                logging.debug(f"Region classified as {region_type}, not intergenic, retrying...")
                continue
            
            # Try to find suitable splice sites
            circ_start, circ_end = gtf_parser.create_suitable_splice_sites(
                genome, chrom, circ_start, circ_end, strand, splice_type)
            
            # Generate the sequence
            circ_seq = genome.get_sequence(chrom, circ_start-1, circ_end, strand)
            
            # Verify if sequence was obtained
            if not circ_seq:
                logging.warning(f"Failed to obtain sequence for intergenic circRNA at {chrom}:{circ_start}-{circ_end}")
                continue
            
            # Generate random copy number for rolling circle structure
            copy_number = self.generate_copy_number()
            
            # Create CircRNA entry - for intergenic, create a synthetic ID with not_in_ENS flag
            # Format: ENSXXXXX|not_in_ENS|intergenic|chrom:start-end|strand|length
            fake_ens_id = f"ENST{np.random.randint(10000, 99999)}0001"
            circ_id = f"{fake_ens_id}|not_in_ENS|intergenic|{chrom}:{circ_start}-{circ_end}|{strand}|{len(circ_seq)}"
            
            # Create pseudo-circular sequence with multiple copies
            pseudo_seq = self.generate_rolling_circle_sequence(circ_seq, copy_number)
            
            # Random expression level (using beta distribution for skewed values)
            tpm = np.random.beta(0.2, 3.0) * 40  # Lowest expression
            
            results.append({
                'id': circ_id,
                'chrom': chrom,
                'start': circ_start,
                'end': circ_end,
                'strand': strand,
                'type': 'intergenic',
                'gene_id': 'NA',
                'tx_id': 'NA',
                'circ_seq': circ_seq,
                'pseudo_seq': pseudo_seq,
                'mature_length': len(circ_seq),
                'exon_count': exon_count,
                'copy_number': copy_number,
                'tpm': tpm,
                'splice_type': splice_type
            })
        
        if attempts >= max_attempts:
            logging.warning(f"Hit maximum attempts ({max_attempts}) while generating intergenic circRNAs")
        
        return results
    
    def generate_rolling_circle_sequence(self, circ_seq, copy_number):
        """Generate a pseudo-circular sequence with multiple copies."""
        # Random starting position in the circRNA
        if len(circ_seq) > 1:
            start_pos = np.random.randint(0, len(circ_seq))
        else:
            start_pos = 0
        
        # Generate the rolling circle sequence
        pseudo_seq = circ_seq[start_pos:] + circ_seq * copy_number
        return pseudo_seq
    
    def simulate_circrnas(self, num_circrnas):
        """Simulate circRNAs based on the specified distribution."""
        # Calculate how many of each type to simulate
        num_eciRNA = int(num_circrnas * self.type_distribution['eciRNA'])
        num_EIciRNA = int(num_circrnas * self.type_distribution['EIciRNA'])
        num_ciRNA = int(num_circrnas * self.type_distribution['ciRNA'])
        num_intergenic = int(num_circrnas * self.type_distribution['intergenic'])
        
        # Adjust to ensure total is correct
        total = num_eciRNA + num_EIciRNA + num_ciRNA + num_intergenic
        if total < num_circrnas:
            num_eciRNA += (num_circrnas - total)
        
        logging.info(f"Simulating {num_circrnas} circRNAs: {num_eciRNA} eciRNA, {num_EIciRNA} EIciRNA, " +
                    f"{num_ciRNA} ciRNA, {num_intergenic} intergenic")
        
        # Simulate each type
        eciRNA_results = self.simulate_eciRNA(num_eciRNA, self.genome, self.gtf_parser)
        EIciRNA_results = self.simulate_EIciRNA(num_EIciRNA, self.genome, self.gtf_parser)
        ciRNA_results = self.simulate_ciRNA(num_ciRNA, self.genome, self.gtf_parser)
        intergenic_results = self.simulate_intergenic_circRNA(num_intergenic, self.genome, self.gtf_parser)
        
        # Combine all results
        all_results = eciRNA_results + EIciRNA_results + ciRNA_results + intergenic_results
        
        # Log statistics
        logging.info(f"Successfully simulated {len(all_results)} circRNAs")
        logging.info(f"  eciRNA: {len(eciRNA_results)} ({len(eciRNA_results)/len(all_results)*100:.1f}%)")
        logging.info(f"  EIciRNA: {len(EIciRNA_results)} ({len(EIciRNA_results)/len(all_results)*100:.1f}%)")
        logging.info(f"  ciRNA: {len(ciRNA_results)} ({len(ciRNA_results)/len(all_results)*100:.1f}%)")
        logging.info(f"  intergenic: {len(intergenic_results)} ({len(intergenic_results)/len(all_results)*100:.1f}%)")
        
        return all_results
    
    def write_output_files(self, circrnas):
        """Write simulated circRNAs to output files."""
        # Create output directory
        os.makedirs(self.output_dir, exist_ok=True)
        
        # Output FASTA file with pseudo-circular sequences
        fa_file = os.path.join(self.output_dir, "circRNAs.fa")
        logging.info(f"Writing pseudo-circular sequences to {fa_file}")
        
        with open(fa_file, 'w') as f:
            for circ in circrnas:
                f.write(f">{circ['id']}\n")
                # Write in FASTA format (line wrap at 60 characters)
                seq = circ['pseudo_seq']
                for i in range(0, len(seq), 60):
                    f.write(f"{seq[i:i+60]}\n")
        
        # Output abundances file with mature circRNA sequences and TPM values
        abund_file = os.path.join(self.output_dir, "abundances.tsv")
        logging.info(f"Writing mature sequences and abundances to {abund_file}")
        
        with open(abund_file, 'w') as f:
            f.write("circRNA_id\tmature_sequence\tTPM\n")
            for circ in circrnas:
                f.write(f"{circ['id']}\t{circ['circ_seq']}\t{circ['tpm']:.2f}\n")
        
        # Output metadata file with detailed information about each circRNA
        meta_file = os.path.join(self.output_dir, "circRNA_metadata.tsv")
        logging.info(f"Writing metadata to {meta_file}")
        
        with open(meta_file, 'w') as f:
            f.write("circRNA_id\tchrom\tstart\tend\tstrand\ttype\tgene_id\ttx_id\t" +
                   "mature_length\texon_count\tcopy_number\tsplice_type\tTPM\n")
            for circ in circrnas:
                # Ensure Ensembl IDs in the output
                gene_id = circ['gene_id'] if circ['gene_id'].startswith('ENS') else circ['gene_id']
                tx_id = circ['tx_id'] if circ['tx_id'].startswith('ENS') else circ['tx_id']
                
                f.write(f"{circ['id']}\t{circ['chrom']}\t{circ['start']}\t{circ['end']}\t" +
                       f"{circ['strand']}\t{circ['type']}\t{gene_id}\t{tx_id}\t" +
                       f"{circ['mature_length']}\t{circ['exon_count']}\t{circ['copy_number']}\t" +
                       f"{circ['splice_type']}\t{circ['tpm']:.2f}\n")
        
        # Output BED file for visualization
        bed_file = os.path.join(self.output_dir, "circRNAs.bed")
        logging.info(f"Writing BED file to {bed_file}")
        
        with open(bed_file, 'w') as f:
            for circ in circrnas:
                name = circ['id']
                score = int(min(1000, circ['tpm'] * 10))  # Scale TPM to 0-1000 range
                
                # RGB color based on type
                rgb = "0,0,0"  # Default black
                if circ['type'] == 'eciRNA':
                    rgb = "27,158,119"  # Teal
                elif circ['type'] == 'EIciRNA':
                    rgb = "217,95,2"    # Orange
                elif circ['type'] == 'ciRNA':
                    rgb = "117,112,179"  # Purple
                elif circ['type'] == 'intergenic':
                    rgb = "231,41,138"   # Pink
                
                f.write(f"{circ['chrom']}\t{circ['start']-1}\t{circ['end']}\t{name}\t{score}\t{circ['strand']}\t" +
                       f"{circ['start']-1}\t{circ['end']}\t{rgb}\t1\t{circ['mature_length']}\t0\n")
        
        logging.info("Successfully wrote all output files")
        
        return fa_file, abund_file, meta_file, bed_file

def main():
    parser = argparse.ArgumentParser(description='Simulate circRNAs with realistic biological features')
    parser.add_argument('--gtf', required=True, help='GTF file with gene annotations')
    parser.add_argument('--genome', required=True, help='Reference genome in FASTA format')
    parser.add_argument('--output', required=True, help='Output directory')
    parser.add_argument('--count', type=int, default=1000, help='Number of circRNAs to simulate (default: 1000)')
    parser.add_argument('--bedtools', default='bedtools', help='Path to bedtools executable (default: bedtools)')
    args = parser.parse_args()
    
    # Set up logging
    setup_logging(args.output)
    
    try:
        # Load genome
        genome = GenomeLoader(args.genome)
        
        # Parse GTF
        gtf_parser = GTFParser(args.gtf, args.bedtools)
        
        # Initialize simulator
        simulator = CircRNASimulator(gtf_parser, genome, args.output, args.bedtools)
        
        # Simulate circRNAs
        circrnas = simulator.simulate_circrnas(args.count)
        
        # Write output files
        simulator.write_output_files(circrnas)
        
        # Clean up
        gtf_parser.cleanup()
        
        logging.info("circRNA simulation completed successfully")
        
    except Exception as e:
        logging.error(f"Error in circRNA simulation: {e}")
        raise
    
if __name__ == "__main__":
    main()