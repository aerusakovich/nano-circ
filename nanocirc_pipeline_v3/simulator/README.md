# nanocirc-sim (simulator v2)

A fast circRNA simulator updated from original publication. It builds realistic circRNA sequences from a genome and a gene annotation (GTF), then checks that simulated long reads actually cross the back-splice junction (BSJ).

You can run it with default settings, or fit it to a real circRNA database first, so the simulated circRNAs look like real ones from that database.

## What it does

- Builds circRNA sequences from a genome and a GTF file.
- Supports 5 circRNA types: `ecirna`, `eicirna`, `cirna`, `antisense`, `intergenic`. See [circRNA types](#circrna-types).
- Can fit its settings to a real circRNA database, so lengths, exon counts, and splice sites match real data.
- Checks that simulated long reads truly cross the BSJ, not just that they came from a "circular" simulation run.
- Every run with the same `--seed` gives the exact same output. Good for reproducible benchmarks.

## Install

```bash
conda env create -f environment.yml
conda activate nanocirc-sim
pip install -e .
```

Or run it straight from the source tree, no install needed:

```bash
PYTHONPATH=. python -m nanocirc_sim.cli ...
```

## Try it now

```bash
./run_example.sh   # builds a tiny toy genome and runs simulate on it, in seconds
./run_tests.sh      # runs the test suite (89 tests, no real data needed)
```

Both work right after a fresh checkout. No install and no external data needed.

## How to use

### 1. Simulate circRNAs

You need a genome FASTA and a gene annotation GTF.

```bash
nanocirc-sim simulate \
  --genome genome.fa --gtf annotation.gtf \
  --n-ecirna 100 --n-eicirna 300 --n-cirna 50 --n-antisense 10 --n-intergenic 20 \
  --total-reads 200000 --seed 42 \
  --outdir results/
```

This writes `circRNAs.fa`, `ground_truth.tsv`, and other files. See [Outputs](#outputs-in---outdir).

### 2. Two ways to pick circRNA counts

You can ask for an exact number per type, or a total number split by a realistic mix.

| Mode | Flags | What it does | When to use it |
| --- | --- | --- | --- |
| **Exact counts** | `--n-ecirna`, `--n-eicirna`, `--n-cirna`, `--n-antisense`, `--n-intergenic` | You set the count for each type yourself. | You want the same number of each type. Good for testing tools fairly across types. |
| **Fitted total** | `--total N` | Splits `N` across types using the config's `type_fractions`. | You want a dataset that matches how common each type really is. |

Exact counts always win if you set both. A config still controls each type's length, exon count, and splice sites either way, only the counts differ.

```bash
# Exact: 2000 of each type
nanocirc-sim simulate --genome genome.fa --gtf annotation.gtf \
  --n-ecirna 2000 --n-eicirna 2000 --n-cirna 2000 --n-antisense 2000 --n-intergenic 2000 \
  --outdir results_exact/

# Fitted: 10,000 total, split by a real database's type mix
nanocirc-sim simulate --genome genome.fa --gtf annotation.gtf \
  --config configs/fitted.yaml --total 10000 \
  --outdir results_fitted/
```

### 3. Fit a config from a real database (optional)

If you have a real circRNA database as a BED12 file, you can fit a config to it. This sets realistic type fractions, lengths, exon counts, and splice sites.

```bash
nanocirc-sim fit-config \
  --db-bed circbase_circatlas_merged.bed12 \
  --genome data/human/hg38.fa --gtf data/human/gencode.v47.gtf \
  --out configs/fitted.yaml
```

Then pass `--config configs/fitted.yaml` to `simulate`. If you skip this step, `simulate` uses the defaults in `configs/default.yaml`.

Types with too few real examples in the database keep the default config's values instead of being fit from noise.

### 4. Check the reads (mandatory for long reads)

If you simulate long reads from `circRNAs.fa` (for example with NanoSim), you must check that "circular" reads really cross the BSJ. Some reads land inside one copy of the sequence and never touch the join.

```bash
nanocirc-sim check-bsj \
  --fastq simulated_reads.fastq \
  --mature-fasta results/circRNAs.mature.fa \
  --metadata results/circRNA_metadata.tsv \
  --outdir results/bsj_check/
```

This writes `bsj_check_reads.tsv` (per read) and `ground_truth_reconciled.tsv` (per circRNA: how many reads were targeted vs. confirmed).

## circRNA types

- `ecirna`: made of exons only. Can skip an internal exon.
- `eicirna`: exons plus one or more retained introns.
- `cirna`: made of a single intron.
- `antisense`: sits over one or more exons of a gene, but on the opposite strand. Uses the gene's exon coordinates, but flips the strand.
- `intergenic`: has no gene overlap on either strand. Strand is reported as `.`, since there is no real gene to pick a strand from.

## Outputs (in `--outdir`)

- `circRNAs.fa`: circular (rolling-circle) sequences. Use this for a long-read simulator, such as NanoSim's `transcriptome` mode.
- `circRNAs.mature.fa`: one linear copy per circRNA, used by `check-bsj`.
- `abundances.tsv`: TPM per circRNA, for NanoSim's expression profile.
- `circRNA_metadata.tsv`: full per-isoform table. Type, gene/transcript/locus, exon/intron counts, splice type, TPM, and target read count.
- `circRNAs.bed`: BED12 with the real exon/intron structure, not one flat block.
- `ground_truth.tsv`: combines the BED12 fields with type, locus/isoform, TPM, and read counts.
- `run_manifest.yaml`: the config used, the seed, input file details, and how many circRNAs were made per type.

## Multiple isoforms per locus

`--isoform-mean` and `--isoform-max` control how often one locus makes more than one circRNA. Isoforms from the same locus share a `locus_id`, but get their own `isoform_index`, and can differ in exon boundaries.

`same_bsj_isoform_prob` in the config sets how often a new isoform reuses an older isoform's exact back-splice boundary, only changing the internal exon pattern. This is a common real-world pattern: one BSJ, several spliced variants.

## Containers

`containers/nanocirc-sim.def` (Apptainer) and `containers/Dockerfile` build the same image. Build from this directory:

```bash
apptainer build containers/nanocirc-sim.sif containers/nanocirc-sim.def
```

`containers/nanosim.def` builds a NanoSim 3.2.3 container, for the read simulation step.

> [!TIP]
> On a cluster where apptainer is not on the login node's PATH, check a compute node first (`srun -p <partition> --pty apptainer build ...`) before assuming it is missing.

## Simulating reads on a cluster (SLURM)

Three steps, each its own sbatch script. `run_read_simulation.sbatch` and `run_pool_select_and_filter.sbatch` chain automatically. `run_pool_select_and_filter.sbatch` also submits `run_linear_matched.sbatch` once per replicate run:

```
                                                                       --> run_linear_matched.sbatch (run 1)
run_read_simulation.sbatch  -->  run_pool_select_and_filter.sbatch    --> run_linear_matched.sbatch (run 2)
   (one big mega-pool)              (pick reads, draw N replicate     --> run_linear_matched.sbatch (run 3)
                                      runs)                               (one per replicate)
```

Run `scripts/build_transcriptome.py` once per genome first, to build the spliced transcriptome FASTA that NanoSim needs.

### 1. Mega-pool

Simulates a large pool of reads with NanoSim, using a real FASTQ to learn realistic read length and error rates, then checks the BSJ.

```bash
python scripts/build_transcriptome.py --gtf annotation.gtf --genome genome.fa --out transcriptome.fa

sbatch --export=ALL,\
DATASET_DIR=$PWD/datasets/mouse_dataset1,\
GENOME_FA=/path/to/genome.fa,\
GTF=/path/to/annotation.gtf,\
TRANSCRIPTOME_FA=$PWD/datasets/mouse_transcriptome.fa,\
CONTROL_FASTQ=/path/to/real_ont_reads.fq,\
N_READS=5000000,\
OUTDIR=$PWD/datasets/mouse_dataset1/mega,\
NANOSIM_SEED=1 \
  scripts/run_read_simulation.sbatch
```

`CONTROL_FASTQ` must be plain FASTQ (not gzipped), and should come from the same real protocol you are simulating. NanoSim uses it to learn its read length and error model, it does not copy its sequences.

`N_READS` should be well above your target read count (5x is a good default). Not every read confirms the BSJ, so a bigger pool means fewer circRNAs end up short.

### 2. Pick reads and draw replicate runs

Picks up to each circRNA's target read count from the mega-pool, then draws `N_RUNS` replicate datasets from that pool.

```bash
sbatch --export=ALL,\
DATASET_DIR=$PWD/datasets/mouse_dataset1,\
GENOME_FA=/path/to/genome.fa,\
TRANSCRIPTOME_FA=$PWD/datasets/mouse_transcriptome.fa,\
CONTROL_FASTQ=/path/to/real_ont_reads.fq,\
CIRC_OUTDIR=$PWD/datasets/mouse_dataset1/mega,\
OUTDIR=$PWD/datasets/mouse_dataset1/select,\
FINAL_OUTDIR=$PWD/datasets/mouse_dataset1/final_balanced,\
N_RUNS=3,READS_PER_RUN=300000 \
  scripts/run_pool_select_and_filter.sbatch
```

Each replicate run picks a different, mostly-overlapping set of circRNAs, using a Poisson draw around each circRNA's expected share. This mimics real sequencing replicates: high-abundance circRNAs show up in every run, low-abundance ones drop in and out.

Each run writes `run_circrnas.bed` and `run_ground_truth.tsv`, listing only the circRNAs that run actually got BSJ-confirmed reads for (see `check-bsj` above), not just circRNAs a read was simulated for. A circRNA with zero confirmed reads in a run is left out of that run's ground truth entirely, it is not written with a zero count. A detection tool cannot be fairly judged against a circRNA that had no real BSJ evidence in the data it was given.

### 3. Add matched linear background reads

Simulates linear (non-circular) reads from the real transcriptome, one for each circular read in a run. This gives every replicate a 1:1 mix of circular and linear reads to challenge the tools, since real circRNA-enriched protocols still have linear contamination.

`run_pool_select_and_filter.sbatch` submits this step for you, once per replicate.

All three sbatch scripts use `#SBATCH --exclusive`. Two NanoSim jobs on the same node can silently crash each other, so each job gets a full node to itself.

## Tests

```bash
./run_tests.sh
# or, if installed: pytest tests/
```

89 unit tests, no real genome or GTF needed. `tests/conftest.py` builds small synthetic genomes on the fly, and each test checks one behavior (for example, that antisense circRNAs land on the flipped strand). Runs in under 2 seconds.
