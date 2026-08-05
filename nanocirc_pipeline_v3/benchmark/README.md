# nanocirc-benchmark

Benchmark for long-read circRNA callers and quantifiers. It scores structural precision/recall, read-count correlation, and detection sensitivity against ground truth.

It does not run any caller itself. It reads outputs already produced by nf-core/nanocirc and the 4 competing tools (isoCirc, circFL-seq, CIRI-long, CircNick-LRS), and scores them.

## History

This is the third version of the benchmark analysis in this repo that was motivated by development of nanocirc pipeline and closing the limitations of the previous benchmark.

`manuscript_v1/5_performance_metrics/` had 7 separate scripts (one for runtime/memory, one for type and length distributions, one for UpSet overlap plots, one for tool-pair combinations, one for precision/recall/F1, one for expression correlation). Each had its own hardcoded genome, GTF, and BED12 paths, and only compared 3 raw tools (CIRI-long, isoCirc, CircNick) with no config file.

`manuscript_v2` moved simulation, tool execution, and BED12 conversion into a Nextflow pipeline for reproducibility. Resource usage (runtime/memory) moved into a Nextflow module (`performance.nf`), reading Nextflow's own trace file. The actual scoring (precision/recall, type, length, UpSet, expression) stayed a separate shell and R/Python script (`Visualisation/plots.sh`) with the same hardcoded paths as v1, and still only the same 3 raw tools.

Compared to our earlier benchmark (Rusakovich et al., 2026), this version (v3) differs in four ways:

1. We found a real bug in `bedtools intersect` when `--split` and `-f`/`-r` are used together (details below). This bug made the earlier benchmark undercount exon-aware detection for every tool. Please consult v3 plots for renewed benchmark plots.
2. All four tools now run from container images, with upstream bug fixes applied (see Containers below).
3. circFL-seq is scored for the first time, alongside the other three tools.
4. Structural matching now also includes BSJ-based matching, both exact (`strict_bsj`) and with tolerance (`relaxed_bsj`, tested at ±5bp and ±10bp), not just exon-based matching.

### The bedtools `--split`/`-f` bug

`bedtools intersect --split -f X -r` looks correct on paper. It runs with no error or warning and gives plausible-looking output. But once more than one `-b` record can overlap the same `-a` record, which is normal in real annotation files containing isoforms, it has a failure mode.

The bug happens when a query record has more than one candidate match nearby. bedtools first finds every `-b` record whose genomic span touches the `-a` record's span, exon-aware or not. Then, seemingly, instead of checking each candidate on its own, it adds up the block lengths of every candidate into one shared denominator and checks a single combined fraction for the whole group. If that combined number misses the threshold, every pair in the group is dropped, including ones that would pass fine by themselves.

#### To reproduce it

No special data needed, just two small BED12 files and bedtools (tested on 2.30.0):

```
$ cat > a.bed <<'EOF'
chr1    1000    2000    A1    0    +    1000    2000    0    1    1000,    0,
EOF

$ cat > b.bed <<'EOF'
chr1    1000    2000    B1_perfect_match    0    +    1000    2000    0    1    1000,    0,
chr1    1000    6000    B2_different_isoform    0    +    1000    6000    0    2    1000,2000,    0,3000,
EOF
```

A is a single 1000bp exon. B1 is an exact copy of A, same single exon, same coordinates. B2 is a different, longer isoform: it shares that same 1000bp exon with A, but also has an unrelated 2000bp exon further downstream, for 3000bp of exons in total. Neither A nor B1 nor B2 has anything to do with a threshold as low as 0.90, they should just match.

First, intersect A against B1 only:

```
$ grep B1_perfect_match b.bed > b1_only.bed
$ bedtools intersect -a a.bed -b b1_only.bed -split -r -f 0.90 -wa -wb
```

Output: one line. A and B1 match perfectly, fraction 1000/1000 = 1.000, as expected.

But if we run the exact same command against the full `b.bed`, with B2 included:

```
$ bedtools intersect -a a.bed -b b.bed -split -r -f 0.90 -wa -wb
```

Output: empty. B1 disappears too, just because B2 is also in the file.

Why we think it happens: bedtools adds B1's length (1000bp) and B2's length (3000bp) into one shared denominator, 4000bp, and checks a single combined fraction for the whole group: 1000/4000 = 0.25. That fails at f=0.90, so bedtools drops the whole group, including the record that was a perfect match. Changing `-f` on this same file confirms the exact breakpoint: it passes up to `f=0.25` and fails from `f=0.26` on, matching 1000/4000.

To confirm this is the actual mechanism, we read the bedtools2 source (`BlockMgr::findBlockedOverlaps`, in `src/utils/FileRecordTools/Records/BlockMgr.cpp`) and built small test cases where we knew the expected number in advance. Two identical B records at the same spot (829bp each) break the match once the threshold goes above 829/1658 = 0.500, exactly as predicted from the source. Three identical copies break it above 829/2487 = 0.333.

This bug is specific to `-split` use and doesn't reproduce itself with the `-split` dropped. It has been reported as an issue `#1142` in `arq5x/bedtools2`.

On one real comparison (CIRI-long vs CircNick-LRS BED12, one mouse sample, threshold 0.95):

| Metric | Value |
| --- | --- |
| Raw overlapping pairs (before threshold) | 19,905 |
| Correct method (`-wo` + hand calculation), passing pairs | 663 |
| bedtools `-f -r`, A→B passing pairs | 359 |
| bedtools `-f -r`, B→A passing pairs | 424 |
| Agreement between methods | 357 |
| Pairs found in only one direction | 139 |
| False negatives from the bug | 306 (46%) |
| False positives from the bug | 2 |
| Jaccard agreement between methods | 53.7% |

A few of the 306 discarded pairs, all of which are exact or near-exact matches by the correct method, and all of which have more than one candidate sharing their locus:

| A record | B record | correct fraction | bedtools `-r -f 0.95` |
| --- | --- | --- | --- |
| chr10:105369950-105370778 | mmu_circ_1788_Tmtc2 | 1.000 | dropped (3 candidates share this locus) |
| chr11:104350487-104378826 | mmu_circ_0341_Kansl1 | 1.000 | dropped (6 candidates share this locus) |
| chr10:12436281-12455564 | mmu_circ_0565_Utrn | 1.000 | dropped (5 candidates share this locus) |

Symmetry: since the group of candidates that get summed together depends on whichever file is `-b`, swapping `-a` and `-b` changes the result for the same pair:

| A record | B record | A -> B | B -> A |
| --- | --- | --- | --- |
| chr10:6828885-6830497 | mmu_circ_2346_Oprm1 | pass | fail |
| chr10:109764518-109765006 | mmu_circ_5638_Nav3 | fail | pass |
| chr11:106096991-106123906 | mmu_circ_3369_Map3k3 | pass | fail |

All three pairs above have identical coordinates and identical exon structure on both sides. There is no biological reason for direction to matter here.

Two more problems, beyond the false negatives. First, because bedtools decides pass or fail for the whole group of candidates at once, a pair that only partly covers its match can still get reported, as long as some other candidate at the same spot makes the group's combined fraction pass. Second, the result is not symmetric, for the reason in the table above: whichever file is given as `-b` is the one whose records get grouped and summed, so `-a`/`-b` order changes the outcome even for a pair that matches perfectly in both directions.

The fix, what v3 does: run `bedtools intersect --split -wo` to get the raw overlap in base pairs, then compute the fraction ourselves, one pair at a time. This checks every pair on its own, with no grouping, so it is symmetric and always uses each record's own real block length as the denominator.


### Containers

None of the four tools originally came with a container. Docker and Singularity images for all four are built and published at https://github.com/aerusakovich/containers_circRNA_tools (DOI: 10.5281/zenodo.19202962), on Quay.io, tested on Linux x86_64. Several images also fix real upstream bugs that otherwise cause silent failures or crashes:

- `quay.io/anrusakovich/ciri-long`: Python 3.8, CIRI-long installed with pip. Giving it too many CPU cores can cause out-of-memory errors, if that happens, use fewer cores.
- `quay.io/anrusakovich/isocirc`: Python 3.11, bedtools ≥2.27.0, minimap2 ≥2.11.
- `quay.io/anrusakovich/circfl-seq`: Python 3.9, minimap2, samtools, bedtools, TRF, TideHunter, plus 4 bug fixes.
- `quay.io/anrusakovich/circnick-lrs`: Python 3.7.6, bedtools 2.29.2, samtools 1.9, pblat 2.5, NanoFilt, deeptools, plus 2 bug fixes. Ships with pre-built human and mouse reference indexes, so no reference setup is needed at runtime.

CIRI-long and circFL-seq need AVX2-compatible hardware (`--constraint=avx2` in Slurm). On HPC and Slurm, use Singularity, not Docker.

### Other method changes in v3

Beyond the bedtools fix above, the scoring method changed in a few more ways:

- **Strand check.** v1 and v2's precision/recall/F1 scoring never checked strand. A call on the wrong strand still counted as a true positive if the coordinates overlapped enough. v3 requires the strand to match as well, except for intergenic calls, which have no biologically meaningfull strand.
- **antisense type.** v1 and v2's ground truth had 4 circRNA types (eciRNA, EIciRNA, ciRNA, intergenic), no antisense type. Now as simulator expanded to this type, evaluation also has been added. 
- **Matching strategies.** v1 and v2 only had one kind of threshold: an arbitrary reciprocal overlap fraction (25%, 50%, 75%, or 95%, picked per script). v3 adds a boundary-tolerance match (`relaxed_bsj`: both ends independently within N bp of ground truth, tested at ±5bp and ±10bp) and an exact match (`strict_bsj`), on top of the corrected exon-overlap match (`exon_based`) described above.
- **Read-count matching.** To decide which tool call and which ground-truth locus are "the same circRNA" for the expression correlation plots, v1 used a fixed 95% bedtools reciprocal-overlap match, the same kind of check found buggy for TP/FP/FN scoring (see History). v3 instead reuses its own boundary-tolerance match (both ends within N bp of ground truth, the same check `relaxed_bsj` uses) for this too, so there is no overlap-fraction threshold to pick at all.
- **Tools scored.** v1 and v2 only scored CIRI-long, isoCirc, and CircNick. v3 adds circFL-seq and nanocirc pipeline's 4 confidence tiers.

Packaging changed too: one config file instead of hardcoded paths in each script, and its own BED12 I/O and bedtools wrapper built in, as well as tests for all used methods.

## Install

```bash
conda env create -f environment.yml
conda activate nanocirc-benchmark
pip install -e .
```

Or run it straight from the source tree, no install needed:

```bash
PYTHONPATH=. python -m nanocirc_benchmark.cli ...
```

## Quick start

```bash
nanocirc-benchmark all \
  --dataset-root /path/to/benchmark_dataset \
  --species mouse \
  --config configs/default.yaml
```

This runs every metric module, in order, and writes every table and plot under `<dataset_root>/benchmark/`. Add `--gtf path/to/annotation.gtf` to also run type classification. Skip it, and that one module is left out.

Each module can also run on its own:

```bash
nanocirc-benchmark performance --dataset-root ... --species mouse --config ...
```

Subcommands: `performance`, `expression`, `sensitivity`, `type-classification`, `overlap`, `length`, `reproducibility`.

See `scripts/run_benchmark.sbatch` for a ready-to-submit SLURM version of the command above.

## Expected dataset layout

```
<dataset_root>/
├── fastq/<run>.fastq.gz
├── ground_truth/<run>_ground_truth.tsv
├── results/<species>/circrna/<run>/...            per-run tool + NanoCirc outputs
├── results/<species>/circrna/crossrun/<species>/<tier>/<species>_<tier>_crossrun.bed12
└── benchmark/
    ├── gt_prefiltered/<run>_gt_prefiltered.bed12   GT restricted to loci with a confirmed circular read in this run
    └── <species>.genome                            bedtools-style chrom-sizes file
```

`fastq/`, `ground_truth/`, `results/`, `gt_prefiltered/`, and `<species>.genome` must already exist before you run this package. They come from the simulator or your own real-data setup, not from this package. Everything else under `benchmark/` is written by this package.

## Metric modules

- **`performance`**: structural precision/recall/F1 for every caller (4 tools plus NanoCirc's 4 confidence tiers), every run, and 3 match strategies (`strict_bsj`, `relaxed_bsj`, `exon_based`). Writes `benchmark_metrics_<species>_<tag>.csv` and a precision/recall scatter plot.
- **`expression`**: read-count correlation (Pearson, Spearman, R²) between each caller and ground truth, for the same 3 match strategies as `performance`. Writes one CSV and one scatter/heatmap plot per strategy: `expression_correlation_<species>.csv`/`.png` for `relaxed_bsj` (the long-standing default, unsuffixed for backward compatibility), plus `_strict_bsj` and `_exon_based` suffixed versions.
- **`sensitivity`**: does a caller find the circRNA at all, split by ground-truth expression level (Low/Medium/High). Writes a CSV and plot.
- **`type-classification`**: classifies every circRNA into one of 5 types (`ecirna`, `eicirna`, `cirna`, `antisense`, `intergenic`), the same way nanocirc's own pipeline does. Needs `--gtf`. Writes a CSV and plot per match strategy.
- **`overlap`**: an UpSet-style plot showing which circRNAs each caller finds, and where callers agree or disagree. Makes both an exon-based and a relaxed-BSJ version, since a caller can score very differently on each.
- **`length`**: mature (spliced) length of true-positive calls, by caller. Also makes a TP/FP/FN heatmap, reusing `performance`'s own output (run `performance` first).
- **`reproducibility`**: a pairwise scatter matrix comparing read counts across every caller and tier, including ground truth, for the same 3 match strategies as `performance`. Writes `reproducibility_pairs_<species>.png` for `relaxed_bsj` (unsuffixed default), plus `_strict_bsj` and `_exon_based` suffixed versions.

## Notes on how matching works

- **Coordinate matching uses a tolerance, not an exact match.** A caller's own called boundary is rarely identical, base for base, to ground truth. `matching.py` matches within a tolerance everywhere, not just in one place.
- **Near-duplicate calls are collapsed, not summed.** If two catalog entries both land within tolerance of one true locus, they are treated as the same locus. Their read count is the max of the two, not the sum. Summing would double-count one real signal as two.
- **Every raw overlap or distance check goes through bedtools, never a Python reimplementation of interval geometry.** Each of the 3 strategies uses a different bedtools call: `strict_bsj` checks exact reciprocal containment (`-f 1.0 -r`), `relaxed_bsj` slops each ground-truth boundary point by the tolerance and intersects it against the caller's own boundary point, `exon_based` gets raw overlap base pairs (`-wo --split`) and computes the fraction itself in Python. Only `exon_based` ever needed that last step, see the bedtools `--split`/`-f` bug in History for why calling bedtools' own `-f`/`-r` fraction flags directly was not safe for that one.
- **Only true-positive calls are used for "what does a call look like" questions**, such as length or type. A call's length only means something once you know it matches a real circRNA.
- **CircNick-LRS read counts are joined through its own BED12 file, not its native output file.** On a run needing `--circnick_liftover_chain`, the native file's own coordinates are still pre-liftover, so they do not match anything else. The BED12 file already has the corrected, post-liftover coordinates.
- **Ground truth only contains circRNAs with a real confirmed circular read in that run.** `gt_prefiltered` is already restricted upstream (in the simulator, at the pool-select stage) to loci with at least one read that was checked and confirmed to actually cross the BSJ, not just a read that was simulated and labeled circular. A circRNA with zero confirmed reads in a run is left out of ground truth entirely, so no tool is scored against a locus it had no real evidence for in the data it was given.

## FP and FN characterization scripts

`scripts/fp_characterize_step1_match_gt.py` and `scripts/fp_characterize_step2_trace_reads.py` trace every false-positive call back to its supporting reads, and show what those reads really are (a real circular read at the wrong strand, at a different locus, boundary noise, or not circular at all). This is the analysis behind the "What are these false positives?" section of `docs/benchmark.md` (in the `nf_nanocirc_long_read` repo).

`scripts/fn_characterize.py` compares every missed ground-truth circRNA against every found one, on length, copy number, expression, splice site, and gene context. This is the analysis behind the "What discovery still misses?" section.

All three scripts need a real benchmark `dataset_root`, plus the matching `nanocirc_sim` `circRNA_metadata.tsv` and a UCSC RepeatMasker BED file. Run any of them with `--help` to see the full argument list.

## Tests

```bash
pytest tests/
```

Tests cover the config loader, the coordinate-matching and near-duplicate-collapse logic, and path resolution. They do not need a real dataset. Full end-to-end runs need one, see `scripts/run_benchmark.sbatch`.
