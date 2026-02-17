# Benchmarking circRNA Detection Tools from Long-Read Sequencing Using Data-Driven and Flexible Simulation Framework

![Version](https://img.shields.io/badge/version-1.0.0-blue)
![License](https://img.shields.io/badge/license-GPL--3.0-green)


## Overview
This is a comprehensive bioinformatics pipeline for generating realistic simulated circular RNA (circRNA) datasets for Oxford Nanopore long-read sequencing. This framework integrates molecular features of circRNAs extracted from established databases and real datasets into NanoSim to create FASTQ reads that accurately reflect biological diversity and technical properties of circular RNAs. CircSim serves as a standardized evaluation framework for benchmarking circRNA detection tools from long-read sequencing data.

![Computational Workflow for Nanopore Sequencing circRNA Simulation](https://drive.google.com/uc?export=view&id=1GEqbbmgli_PldVW9ClEypWZ2WXqb1mSj)
The workflow integrates wet-lab Oxford Nanopore Technologies (ONT) read characteristics with computational simulation techniques. Panels A-C depict sequential stages of data processing and simulation: (A) initial feature extraction from circRNA databases and experimental reads, (B) generation of protocol-specific reference circRNA sequences, and (C) NanoSim-based simulation generating detailed computational read models. 

## Features

- **circRNA simulation**: Generates in silico reads that mimic the molecular and sequencing characteristics of circular RNAs based on real Nanopore sequencing data
- **circRNA types generation**: Simulates four distinct circRNA types (exonic circRNAs, circular intronic RNAs, exon-intron circRNAs, and intergenic circRNAs) with realistic distributions
- **Biological features simulation**: Incorporates key biological features including splice site motifs, exon count distributions, and mature length variations extracted from circRNA databases
- **Customizable parameters**: Offers user-defined controls for rolling circle amplification, read count, and other simulation parameters to match specific experimental conditions
- **Standardized benchmark framework**: Provides a controlled environment with precisely known circRNA annotations for  assessment of detection tool performance
- **Performance evaluation**: Includes descriptive and performance metrics (precision, recall, F1 score) assessment across different overlap thresholds and circRNA subtypes

# Repository Structure

```
nano-circ/
├── envs/                    # Conda environment definitions
├── manuscript_v1/           # Scripts used in manuscript version 1 (shell/Python-based)
│   ├── 0_file_processing/
│   ├── 1_feature_extraction/
│   ├── 2_simulation_circRNAs/
│   ├── 3_tools/
│   ├── 4_make_bed12/
│   └── 5_performance_metrics/
└── manuscript_v2/           # Scripts used in manuscript version 2 (Nextflow pipeline)
    └── Scripts/
        ├── Analysis/        # Nextflow pipeline (main.nf + modules)
        └── Visualisation/   # Updated visualisation scripts
```

> **Note:** The `manuscript_v2` analysis is replacing Nanosim simulation step (step 2.3), tool execution(step 3) and bed12 conversion (step 4) from v1 with a Nextflow pipeline for improved reproducibility and scalability across multiple datasets. The `manuscript_v2` visualisation replaces previous visualisations (step 5) by interacting with nextflow files to calculate performance metrics and resource usage across several runs.

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

## Manuscript v2 Workflow (Nextflow Pipeline)

Manuscript v2 introduces a fully integrated [Nextflow](https://www.nextflow.io/) pipeline that automates the complete analysis — from NanoSim read simulation through BED12 generation and performance evaluation — along with updated visualisation scripts.

### Running the Pipeline

```bash
cd manuscript_v2/Scripts/Analysis
nextflow run main.nf -params-file params.yaml -c nextflow.config
```

Pipeline parameters (reference genome, annotation, read counts, tool paths, etc.) are configured in `params.yaml`. Resource allocation and executor settings are defined in `nextflow.config`.

### Pipeline Modules

The pipeline (`main.nf`) is composed of the following modules located in `manuscript_v2/Scripts/Analysis/modules/`:

| Module | Description |
|---|---|
| `nanosim.nf` | Characterises real ONT reads and simulates circRNA long reads with realistic error profiles |
| `cirilong.nf` | Runs CIRI-long detection and collapses redundant circRNA calls |
| `isocirc.nf` | Runs IsoCirc detection for full-length circRNA isoform characterisation |
| `circnick.nf` | Runs circNICK-lrs detection using BLAT-based alignment |
| `bed12_conversion.nf` | Converts tool-specific outputs and ground truth annotations to standardised BED12 format |
| `performance.nf` | Computes and plots resource usage and execution time based on trace Nextflow files |

### Visualisation

Updated visualisation scripts are located in `manuscript_v2/Scripts/Visualisation/`:

```bash
# Extract required files into split run folders using symlinks
bash manuscript_v2/Scripts/Visualisation/symlinks_nf.sh

# Generate all figures used in the manuscript
bash manuscript_v2/Scripts/Visualisation/plots.sh
```

## License

This project is freely available under the GNU General Public License v3.0 (GPL-3.0). See the [GNU GENERAL PUBLIC LICENSE](https://www.gnu.org/licenses/gpl-3.0.en.html) for details.

---

Developed by the IGDR (Institut de Génétique et Développement de Rennes) - UMR 6290 CNRS, University of Rennes