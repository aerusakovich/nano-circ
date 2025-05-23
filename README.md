# Benchmarking circRNA Detection Tools from Long-Read Sequencing Using Data-Driven and Flexible Simulation Framework

![Version](https://img.shields.io/badge/version-1.0.0-blue)
![License](https://img.shields.io/badge/license-GPL--3.0-green)


## Overview
This is a comprehensive bioinformatics pipeline for generating realistic simulated circular RNA (circRNA) datasets for Oxford Nanopore long-read sequencing. This framework integrates molecular features of circRNAs extracted from established databases and real datasets into NanoSim to create FASTQ reads that accurately reflect biological diversity and technical properties of circular RNAs. CircSim serves as a standardized evaluation framework for benchmarking circRNA detection tools from long-read sequencing data.

![Computational Workflow for Nanopore Sequencing circRNA Simulation](https://drive.google.com/uc?export=view&id=1GEqbbmgli_PldVW9ClEypWZ2WXqb1mSj)
The workflow integrates wet-lab Oxford Nanopore Technologies (ONT) read characteristics with computational simulation techniques. Panels A-C depict sequential stages of data processing and simulation: (A) initial feature extraction from circRNA databases and experimental reads, (B) generation of protocol-specific reference circRNA sequences, and (C) NanoSim-based simulation generating detailed computational read models. 

### CircRNA Types and Protocol-Specific Sequence Patterns

![Types of circRNAs](https://drive.google.com/uc?export=view&id=1zTN8EqZ-gAqHNGAPPAdEq4Gz9jXBywTd)

**Supplementary Figure 1. Types of circRNAs.** (A) EcircRNAs contain only exon sequences. (B) EIciRNAs have both exon and intron sequences. (C) ciRNAs contain only intron sequences. (D) Intergenic circRNAs are derived from the sequence between two different protein-coding genes. Adapted from: Zeng, Le & Liu, Longzhou & Ni, Wen-Juan & Xie, Fuhua & Xiaomin, Leng. (2023). Circular RNAs in osteosarcoma: An update of recent studies (Review). International Journal of Oncology. 63. 10.3892/ijo.2023.5571.

![Expected sequence based on wet-lab approach](https://drive.google.com/uc?export=view&id=1DLjgEFbK0kohZGHeE-xtatbMrRaLbIZj)

**Supplementary Figure 2. Expected sequence based on wet-lab approach.** (A) Nicking approach utilised by circNICK-lrs tool provides one circRNA - one read pattern. (B) Rolling circle approach utilised by CIRI-long and isoCIRC have several repeats of circRNA per read creating during rolling circle reverse transcription or rolling circle amplification step.

## Features

- **Realistic circRNA Simulation**: Generates in silico reads that mimic the molecular and sequencing characteristics of circular RNAs based on real Nanopore sequencing data
- **Comprehensive circRNA Types**: Simulates four distinct circRNA types (exonic circRNAs, circular intronic RNAs, exon-intron circRNAs, and intergenic circRNAs) with realistic distributions
- **Biological Accuracy**: Incorporates key biological features including splice site motifs, exon count distributions, and mature length variations extracted from circRNA databases
- **Customizable Parameters**: Offers user-defined controls for rolling circle amplification, read count, and other simulation parameters to match specific experimental conditions
- **Standardized Benchmark Framework**: Provides a controlled environment with precisely known circRNA annotations for rigorous assessment of detection tool performance
- **Performance Evaluation**: Includes comprehensive metrics (precision, recall, F1 score) assessment across different overlap thresholds and circRNA subtypes

## Installation

### Setting Up Conda Environments

The workflow requires several Conda environments to run different components of the pipeline. The main environment (`bed12`) is used as the default, with additional specialized environments for specific tools.

```bash
# Clone this repository and then cd to it
git clone https://gitlab.com/bioinfog/circall/circ_jobim.git
cd circ_jobim

# Create the main environment
conda env create -f envs/main_environment.yml
conda activate bed12

# Create tool-specific environments (as needed)
conda env create -f envs/nanosim_environment.yml
conda env create -f envs/circnick_environment.yml
conda env create -f envs/ciri-long_environment.yml
conda env create -f envs/isocirc_environment.yml
```

### Environment Details

1. **Main Environment (bed12)**: Default environment for running the core pipeline and analysis tools
2. **NanoSim Environment**: For running the NanoSim read simulator
3. **Tool-specific Environments**:
   - circnick_environment
   - ciri-long_environment
   - isocirc_environment

To ensure reproducibility, each tool runs in its own isolated environment with appropriate dependencies.

## Workflow Steps

### 1. Genome Coordinate Conversion (liftover_script.py)

This Python script converts circRNA coordinates from mm9 to mm10 genome assembly, which is necessary for compatibility with current reference genomes and tools.

```bash
# Usage
python liftover_script.py --fasta input_circrnas.fasta --output liftover_output
```

**Features:**
- Automatically downloads the mm9ToMm10 chain file from UCSC
- Extracts coordinates from circRNA FASTA headers using flexible pattern matching
- Performs coordinate liftOver from mm9 to mm10
- Creates an updated FASTA file with mm10 coordinates

### 2. Feature Extraction and Analysis

This component analyzes circRNA databases to extract biological features:

```bash
# Usage
sbatch 1_feature_extraction_db.sh
```

**Features:**
- Processes and intersects circRNA databases (CircAtlas, CircBase)
- Analyzes circRNAs for splice sites and classification
- Generates comprehensive metrics for each database

### 3. Database Visualization

Creates comprehensive plots for circRNA data analysis:

```bash
# Usage
python 2_database_vis.py --csv circRNA_results.csv --output visualization_output
```

**Features:**
- Optimized pie charts for circRNA type and splice site distributions
- Boxplots showing mature length distribution by circRNA type
- Stacked bar charts for exon count distribution

### 4. Wet Lab Feature Extraction

Extracts key parameters from real sequencing data for use in simulation:

```bash
# Usage
python 4_wet_features.py --bed_file mapped_reads.bed --trf_file tandem_repeats.txt --output_dir feature_extraction_results
```

**Features:**
- Analyzes read length distributions from real data
- Extracts rolling circle features (period sizes, copy numbers)
- Identifies tandem repeat characteristics
- Generates parameters for realistic simulation

### 5. Parallel CircRNA Simulation

Runs multiple parallel simulations to generate circRNA sequences:

```bash
# Usage
sbatch 1_10x_sim.sh
```

**Features:**
- Manages simulation of circRNA datasets
- Automatically generates comprehensive summary reports

### 6. CircRNA Sequence Generation

This Python script simulates circRNAs with biologically accurate features:

```bash
# Usage
python 2_circ_fa_generation_exoncount.py --gtf gencode.annotation.gtf --genome reference.fa --output output_dir --count 10000 --bedtools bedtools_path
```

**Features:**
- Generates four distinct circRNA types (eciRNA, EIciRNA, ciRNA, intergenic)
- Produces realistic mature length distributions based on database analysis
- Creates biologically accurate splice site patterns
- Includes rolling circle amplification with realistic copy numbers

### 7. NanoSim Read Simulation

Generates simulated Nanopore reads from the circRNA sequences:

```bash
# Usage
sbatch 3_nanosim_jobim.sh
```

**Features:**
- Performs read characterization from real Nanopore data
- Simulates realistic circRNA long reads with Oxford Nanopore error profiles
- Generates linear mRNA reads as negative controls
- Creates data with realistic expression quantification

### 8. Data Preparation for Tool Evaluation

Prepares the simulated reads for circRNA detection tool evaluation:

```bash
# Usage
sbatch 0_pooling.sh
```

**Features:**
- Combines circRNA and linear RNA reads into a single FASTQ file
- Adds identifiers to read headers to track their origin
- Creates both uncompressed and compressed versions of the data

### 9. CircRNA Detection Tool Running

#### CIRI-long Detection

```bash
# Usage
sbatch 1_ciri-long.sh
```

**Features:**
- Processes Nanopore long reads using the CIRI-long algorithm
- Identifies candidate circular RNAs based on rolling circle signatures
- Collapses similar circRNAs to reduce redundancy
- Integrates with reference annotations

#### IsoCirc Detection

```bash
# Usage
sbatch 2_isocirc.sh
```

**Features:**
- Processes Nanopore long reads using the IsoCirc pipeline
- Specialized in characterizing full-length circRNA isoforms
- Detects and quantifies alternative splicing events in circRNAs

#### CircNick-lrs Detection

```bash
# Usage
sbatch 3_circnick.sh
```

**Features:**
- Processes Nanopore long reads using the circNICK-lrs workflow
- Specifically designed for long-read sequencing data
- Uses BLAT for more sensitive alignment of circRNAs
- Identifies backsplice junctions using specialized algorithms

### 10. BED12 Format Conversion

These scripts convert tool outputs to standardized BED12 format for evaluation:

```bash
# Usage
sbatch 1_ground_truth.sh
sbatch 2_jobim_circnick_bed.sh
sbatch 3_jobim_cirilong_bed.sh
```

**Features:**
- Converts tool-specific outputs to standardized BED12 format
- Maintains all essential circRNA properties in BED12 format
- Enables accurate overlap calculation in evaluation

**Note:** IsoCirc already outputs results in BED12 format (`isocirc_output/isocirc.bed`), so no conversion is required.

### 11. Performance Metrics Analysis

A suite of scripts for comprehensive evaluation of circRNA detection tool performance:

```bash
# Examples
python 1_job_plots.py --time-csv runtime.csv --memory-csv memory.csv --output-dir metrics_output
python 2_type_length_plots.py --ground-truth ground_truth.bed --circrna-db circRNAs.bed --output-dir type_length_analysis
sbatch 3_upset.sh
sbatch 4_metrics_performance.sh
sbatch 5_combo_performance.sh
```

**Features:**
- Computation performance evaluation with runtime and memory usage visualizations
- CircRNA type and length distribution analysis across tools
- UpSet plots to visualize set intersections between tools
- Precision, recall, and F1 score comparisons with and without exon-aware detection
- Tool combination analysis to determine optimal detection approaches

## License

This project is freely available under the GNU General Public License v3.0 (GPL-3.0). See the [GNU GENERAL PUBLIC LICENSE](https://www.gnu.org/licenses/gpl-3.0.en.html) for details.

---

Developed by the IGDR (Institut de Génétique et Développement de Rennes) - UMR 6290 CNRS, University of Rennes