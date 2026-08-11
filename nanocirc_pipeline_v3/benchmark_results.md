# nanocirc: Benchmark Results and Recommendations

This benchmark is adapted from our own earlier benchmark study, Rusakovich et al., "Benchmarking circRNA detection tools from long-read sequencing using a data-driven and flexible simulation framework" (Peer Community Journal, 2026; see [CITATIONS.md](https://github.com/aerusakovich/nf_nanocirc_long_read/blob/main/CITATIONS.md) for the full reference), modified and extended based on our newer findings. 

This page summarizes the benchmark behind nanocirc's four confidence tiers and gives practical guidance on which one to use. The new simulator and benchmark scripts is in this repo, [nano-circ](https://gitlab.com/bingo-igdr/nano-circ), as v3. Results below are from simulated human and mouse datasets, 3 runs each.

Several benchmarks are shown below:

- **[Same circRNAs, three wet-lab protocols](#same-circrnas-three-wet-lab-protocols)**: NanoSim's read-error model is characterized from real CIRI-long/isoCIRC and circFL-seq protocols ONT control reads; circRNA type composition and length/exon/splice distributions are the same across all 3 wet-lab methods, circRNA generation based on circRNAs from a human circbase+circatlas catalog. 
- **[Species benchmark: human vs mouse](#species-benchmark-human-vs-mouse)**: CIRI-long protocol, each species' own circRNA catalog and genome. Any difference here comes from species biology (genome, catalog composition), not the wet-lab protocol, which is fixed to CIRI-long for both.

## Recommendations

- **Choose circFL-seq wet-lab protocol if you can, followed by CIRI-long, for maximum tool compatibility.** Both work with all 4 tools; circFL-seq protocol gave the highest recall, highest F1, and the highest BSJ-confirmation rate of the 3 protocols we tested (64.1% vs CIRI-long's 54.8%), and its read structure also allowed to better resolve isoform structure. isoCirc-protocol reads are the hardest for the tools to work with, and CIRI-long specifically is almost incompatible with them. See [Overall precision/recall](#overall-precisionrecall) and [Isoform detection](#isoform-detection) under [Same circRNAs, three wet-lab protocols](#same-circrnas-three-wet-lab-protocols).
- **Quantification, Sensitivity/Precision and tool reliability are not connected, but equally important.** All 4 tools needed real bug fixes to run correctly at scale, several silent rather than crashing (wrong fastq/genome argument wiring, a BAM filename mismatch that silently dropped intron-coverage output, a BioPython-version crash that only surfaces on certain read alignments). Separately, a tool being good at finding the right circRNA does not mean its own built-in quantifier counts it well, and a native quantifier can't count a circRNA that tool didn't itself report. A tool-independent, remap-based quantification step (as used in nanocirc-pipeline) avoids both problems: the same counting method applies regardless of which tool/tier sourced a locus, and read counts agree well across tools/tiers and against ground truth (see [Quantification](#quantification)). A shared, cross-sample catalog is also recommended for a valid multi-sample count matrix DESeq2 needs, which quantifying each sample against its own separately-called catalog does not give.
- **For isoform-level differential expression, nanocirc pipeline's `balanced_recall` is a reasonable choice.** In this benchmark, the best single raw tool (circFL-seq) fully resolved 11.5% of human and 15.0% of mouse multi-isoform families on its own; `balanced_recall` reached 14.4%/23.5%, `discovery` 16.3%/26.7% (99% reciprocal-overlap match, see [Isoform detection](#isoform-detection)). `balanced_recall` was close to `discovery`'s isoform resolution with higher precision, which is why we lean toward it for isoform-level work; `discovery` remains a reasonable choice for locus-level differential expression, where the exact internal exon structure matters less. These numbers come from our simulation of multi-isoform families (17-161 per dataset), so treat them as directional rather than definitive.
- **Exception: on isoCirc-protocol data, `discovery` tier looked better than `balanced_recall` for isoform-level work in this benchmark.** CircNick-LRS's own structure calls were unreliable in every dataset we tested, and CIRI-long is close to incompatible with isoCirc-protocol reads, so on this protocol only circFL-seq and isoCirc remain as sources of agreement, and they rarely called the exact same structure independently. `balanced_recall`'s isoform-family resolution was 0% in our isoCirc-protocol runs, `discovery` reached about 12%. Of the two raw tools still usable on this protocol, circFL-seq resolved more families alone than isoCirc (2.0% vs 0.0%), though both numbers are low. See [Isoform detection](#isoform-detection).
- **`balanced_precision`/`high_confidence` only contain calls with cross-tool agreement on both BSJ and structure.** circFL-seq's own raw call was the most structurally accurate single tool in every dataset in this benchmark, though this benchmark covers a limited set of species and protocols. See [Isoform detection](#isoform-detection).

## Same circRNAs, three wet-lab protocols

We used one circRNA catalog (same circbase/circatlas composition, human genome, same abundances) and only changed the wet-lab protocol NanoSim is characterized from. Any difference in the results below comes from the wet-lab protocol characteristics.

### Overall precision/recall

<img src="images/benchmark/pr_scatter_equal_ciri_long.png" width="100%"/>
<img src="images/benchmark/pr_scatter_equal_circfl.png" width="100%"/>
<img src="images/benchmark/pr_scatter_equal_isocirc.png" width="100%"/>

Same circRNA catalog, same tools, same tiers, only the protocol the reads were simulated from changes between these 3 plots.

Mean precision (P) / recall (R) / F1 across 3 runs, exon-based matching, per nanocirc tier, plus each individual tool's own raw calls for comparison:


| Tier / tool | CIRI-long protocol | circFL-seq protocol | isoCirc protocol (CIRI-long on) |
| --- | --- | --- | --- |
| IsoCirc alone | P 0.89 / R 0.14 / F1 0.25 | P 0.90 / R 0.19 / F1 0.32 | P 0.87 / R 0.08 / F1 0.15 |
| CircFL-seq alone | P 0.79 / R 0.21 / F1 0.33 | P 0.77 / R 0.25 / F1 0.37 | P 0.74 / R 0.23 / F1 0.36 |
| CIRI-long alone | P 0.74 / R 0.34 / F1 0.46 | P 0.77 / R 0.41 / F1 0.53 | P 0.77 / R 0.03 / F1 0.05 |
| CircNick-LRS alone | P 0.24 / R 0.09 / F1 0.13 | P 0.24 / R 0.10 / F1 0.14 | P 0.22 / R 0.10 / F1 0.14 |
| `discovery` | P 0.58 / R 0.45 / F1 0.50 | P 0.62 / R 0.53 / F1 0.57 | P 0.54 / R 0.28 / F1 0.37 |
| `balanced_recall` | P 0.74 / R 0.42 / F1 0.53 | P 0.77 / R 0.51 / F1 0.61 | P 0.81 / R 0.25 / F1 0.38 |
| `balanced_precision` | P 0.85 / R 0.22 / F1 0.35 | P 0.85 / R 0.28 / F1 0.42 | P 0.85 / R 0.09 / F1 0.17 |
| `high_confidence` | P 0.89 / R 0.10 / F1 0.18 | P 0.91 / R 0.13 / F1 0.23 | P 0.92 / R 0.02 / F1 0.04 |

### Ground truth coverage

<img src="images/benchmark/equal_gt_coverage_stacked_ciri_long.png" width="100%"/>
<img src="images/benchmark/equal_gt_coverage_stacked_circfl.png" width="100%"/>
<img src="images/benchmark/equal_gt_coverage_stacked_isocirc.png" width="100%"/>

Green+grey sum to 100% of ground truth (a caller's recall split), red stacks separately as false predictions.

### circRNA type classification

<img src="images/benchmark/equal_type_classification_ciri_long_exon.png" width="100%"/>
<img src="images/benchmark/equal_type_classification_ciri_long_relaxed_bsj.png" width="100%"/>
<img src="images/benchmark/equal_type_classification_ciri_long_strict_bsj.png" width="100%"/>
<img src="images/benchmark/equal_type_classification_circfl_exon.png" width="100%"/>
<img src="images/benchmark/equal_type_classification_circfl_relaxed_bsj.png" width="100%"/>
<img src="images/benchmark/equal_type_classification_circfl_strict_bsj.png" width="100%"/>
<img src="images/benchmark/equal_type_classification_isocirc_exon.png" width="100%"/>
<img src="images/benchmark/equal_type_classification_isocirc_relaxed_bsj.png" width="100%"/>
<img src="images/benchmark/equal_type_classification_isocirc_strict_bsj.png" width="100%"/>

Same TP/FP-vs-GT breakdown as the GT coverage plot above, split by circRNA biotype.

### Isoform detection

Long-read protocols are used mainly to resolve internal exon structure, not just BSJ position. Among ground-truth loci where the same BSJ has more than one real isoform (isoform family), we checked how often each caller recovers the exact isoform, versus calls the right BSJ with the wrong structure, versus misses the BSJ entirely.

<img src="images/benchmark/equal_isoform_detection_reciprocal_overlap_99_ciri_long.png" width="100%"/>
<img src="images/benchmark/equal_isoform_detection_reciprocal_overlap_99_circfl.png" width="100%"/>
<img src="images/benchmark/equal_isoform_detection_reciprocal_overlap_99_isocirc.png" width="100%"/>

<img src="images/benchmark/equal_isoform_detection_reciprocal_overlap_90_ciri_long.png" width="100%"/>
<img src="images/benchmark/equal_isoform_detection_reciprocal_overlap_90_circfl.png" width="100%"/>
<img src="images/benchmark/equal_isoform_detection_reciprocal_overlap_90_isocirc.png" width="100%"/>

Green+grey sum to 100% of the isoform-family ground truth; red shows the calls that found the right BSJ but didn't match internal structure criteria. 99% reciprocal exon overlap allows only trivial boundary noise; 90% allows a small sequence mismatch while still requiring the same structure, not a different one.

<img src="images/benchmark/equal_isoform_family_completeness_reciprocal_overlap_99_ciri_long.png" width="100%"/>
<img src="images/benchmark/equal_isoform_family_completeness_reciprocal_overlap_99_circfl.png" width="100%"/>
<img src="images/benchmark/equal_isoform_family_completeness_reciprocal_overlap_99_isocirc.png" width="100%"/>

<img src="images/benchmark/equal_isoform_family_completeness_reciprocal_overlap_90_ciri_long.png" width="100%"/>
<img src="images/benchmark/equal_isoform_family_completeness_reciprocal_overlap_90_circfl.png" width="100%"/>
<img src="images/benchmark/equal_isoform_family_completeness_reciprocal_overlap_90_isocirc.png" width="100%"/>

Per family plot instead of per isoform: did the caller get every isoform in the family right, some of them, or none. Getting one isoform right while missing or misclassifying other is the most common outcome; fully resolving isoform family is a rare event for all tools.

### Expression quantification and reproducibility

<img src="images/benchmark/equal_expression_correlation_ciri_long.png" width="100%"/>
<img src="images/benchmark/equal_expression_correlation_circfl.png" width="100%"/>
<img src="images/benchmark/equal_expression_correlation_isocirc.png" width="100%"/>

Pearson r, Spearman rho, R², read-count correlation against ground truth and cross-method agreement, per tool/tier.

<img src="images/benchmark/equal_sensitivity_by_expression_ciri_long.png" width="100%"/>
<img src="images/benchmark/equal_sensitivity_by_expression_circfl.png" width="100%"/>
<img src="images/benchmark/equal_sensitivity_by_expression_isocirc.png" width="100%"/>

Detection sensitivity by ground-truth expression tertile.

<img src="images/benchmark/equal_reproducibility_pairs_ciri_long.png" width="100%"/>
<img src="images/benchmark/equal_reproducibility_pairs_circfl.png" width="100%"/>
<img src="images/benchmark/equal_reproducibility_pairs_isocirc.png" width="100%"/>

Pairwise read-count agreement across every tool/tier combination, ground truth included.

**Protocol results:**

- **circFL-seq protocol generates reads that work well with all tools**, at every tier: highest recall, highest F1, and the highest BSJ-confirmation rate of the 3 protocols (64.1%).
- **CIRI-long protocol reads are also compatible with all tools**, with a lower BSJ-confirmation rate (54.8%) and lower recall at every tier than circFL-seq.
- **isoCirc protocol reads are the hardest for the tools to work with**, despite carrying the exact same ~9,800 underlying circRNAs. Its recall is the lowest of the 3 protocols at every tier, but precision remains comparable or higher to other wet-lab approaches. CIRI-long is almost incompatible with this wet lab protocol, and that was verified by running CIRI-long on real isoCIRC dataset as well - the tool never managed to finish the analysis on a full dataset.

Because the circRNA catalog, its abundances, and the simulated read amount were identical across all 3 comparisons, this difference comes entirely from how each protocol's wet-lab induced error profile and read length characteristics change the reads NanoSim generates from them.

**Recommendations**: if your protocol choice is flexible, circFL-seq or CIRI-long library prep will give nanocirc's `balanced`/`high_confidence` tiers meaningfully more usable recall than isoCirc protocol. If you already generated your reads and are working with isoCirc data, you can prioritise balanced_recall mode which reaches high precision on this read type and expect potential incopatibility with CIRI-long tool.

circNICK-lrs protocol reads were not tested, as all 3 other tools expect rolling-circle generated data with multiple copies of circRNA per read, making circNICK-lrs protocol compatible only with circNICK-lrs tool, thus it wouldn't benefit from nanocirc's consensus voting algorithm. You still can run circNICK-lrs inside the nanocirc pipeline and benefit from nanocirc's upstream patches to circNICK-lrs tool and nanocirc quantification approach.

CIRI-long was the weakest individual tool on the isoCirc protocol (F1 0.05, far behind its performance on the CIRI-long/circFL-seq protocols above); for a benchmark of the isoCirc protocol with CIRI-long excluded from the tool set entirely (used as nanocirc default for isoCirc wet-lab), see the [nf_nanocirc_long_read pipeline's own benchmark page](https://github.com/aerusakovich/nf_nanocirc_long_read/blob/main/docs/benchmark.md#same-circrnas-three-wet-lab-protocols).

## Species benchmark: human vs mouse

Wet-lab protocol is fixed (CIRI-long) for both species; each species uses its own genome, GTF, and circRNA catalog. Any difference below comes from species circRNA biology and genome/annotation differences, not the wet-lab protocol.

### Overall precision/recall

<img src="images/benchmark/pr_scatter_human.png" width="100%"/>
<img src="images/benchmark/pr_scatter_mouse.png" width="100%"/>

Mean precision (P) / recall (R) / F1 across 3 runs, exon-based matching, per nanocirc tier:

| Tier | Human | Mouse |
| --- | --- | --- |
| `discovery` | P 0.58 / R 0.39 / F1 0.47 | P 0.65 / R 0.54 / F1 0.59 |
| `balanced_recall` | P 0.74 / R 0.36 / F1 0.48 | P 0.85 / R 0.50 / F1 0.63 |
| `balanced_precision` | P 0.85 / R 0.20 / F1 0.32 | P 0.94 / R 0.34 / F1 0.50 |
| `high_confidence` | P 0.91 / R 0.09 / F1 0.17 | P 0.99 / R 0.18 / F1 0.30 |

**What we see:**

- **Mouse outperforms human at every tier**, on both precision and recall. The gap is largest at `balanced_precision` (F1 0.50 vs 0.32) and smallest at `discovery` (F1 0.59 vs 0.47).
- Detailed CIRI-long-protocol breakdowns (ground truth coverage, false-positive characterization, type classification, quantification) for both species are in the [CIRI-long wet-lab protocol benchmark](#ciri-long-wet-lab-protocol-benchmark) section below.

**Recommendations**: species-level differences are real and can be as large as the wet-lab-protocol differences above. Treat published tier-level precision/recall numbers as indicative, not as guarantees that transfer across species, and re-benchmark on your own species/catalog using nanocirc-benchmark if you work with a different one.

## Which tier should I use?

Each tier trades recall against precision differently. There is no single best tier, only the right one for your question.

| Tier | What it optimizes for | Trade-off |
| ---- | ---------------------- | --------- |
| **`discovery`** | Recall. Keeps almost every circRNA any tool found. | Highest false-positive rate: computational artefacts and incorrectly predicted structures will pass alongside real calls. |
| **`high_confidence`** | Precision. Only calls with strong cross-tool agreement on both BSJ and structure or isoCirc calls that met defined read threshold survive. | Sacrifices recall to get the structure right: false positives are minimal, but the reduced set will miss real biology so it's not fit for exploratory studies. |
| **`balanced_recall`** | A middle ground leaning recall. Cuts most of `discovery`'s false-positive rate while keeping reasonable recall. | Good default for exploratory analyses where getting the exact internal exon structure right, or an occasional false positive, isn't critical to the biological question. |
| **`balanced_precision`** | A middle ground leaning precision. Keeps recall relatively high while pushing precision further than `balanced_recall`. | Best for precision-leaning analyses that need to avoid false positives but still want decent recall, not just the aggressively small `high_confidence` set. |

> **Recommendation:** use one of the two **balanced** tiers for most applications. Reserve `discovery` for exploratory, locus-level studies, and `high_confidence` for when you specifically need a small, near-certain set.

## Quantification

Quantification is its own remap step, independent of which tool originally called a locus: the reproducibility and expression correlation plots below show read counts agree well across tools and tiers, and against ground truth, so switching tier for a re-analysis doesn't mean changing to a different counting method.

For differential expression (read counts, DESeq2), `discovery`'s false positives are not easy-to-spot near-zero-read artefacts that low-count filtering could remove (median ~16 reads vs ~16-23 for true positives in our benchmark), so they will be present in the count matrix as real-looking features. `high_confidence`, on the other end, drops enough real calls to hurt power on lower-expressed circRNAs. The Balanced modes avoid both failure modes.

## CIRI-long wet-lab protocol benchmark

Dataset details: NanoSim characterized from real CIRI-long-protocol ONT control reads; circRNA composition fit from a general circbase/circatlas catalog. FP/FN characterization below is human only.

### Overall precision/recall

<img src="images/benchmark/pr_scatter_human.png" width="100%"/>
<img src="images/benchmark/pr_scatter_mouse.png" width="100%"/>

Precision/recall for every tool and every nanocirc tier, at three structural matching strictness levels (exon-based, relaxed BSJ, strict BSJ). The four tiers trace out the trade-off described above: `discovery` sits at high recall/lower precision, `high_confidence` at the opposite corner, and the two balanced tiers in between.

### Ground truth coverage

<img src="images/benchmark/gt_coverage_stacked_human.png" width="100%"/>
<img src="images/benchmark/gt_coverage_stacked_mouse.png" width="100%"/>

Green (recovered) plus grey (missed) always sum to exactly 100% of ground truth, so that split directly shows a caller's recall. Red stacks separately on top of the 100% line as false calls the caller made scaled the same way. 

### What are these false positives?

We traced every FP-supporting read back to the simulator's own ground truth (single-tool discovery exon-based FPs, all 3 human runs, mean % ± SD across runs):

| FP read is... | CIRI-long | CircFL-seq | IsoCirc | CircNick-LRS |
| --- | --- | --- | --- | --- |
| Real circ read, correct boundary, wrong strand (antisense called as sense) | 86.0 ± 0.6 | 79.4 ± 1.7 | 18.8 ± 7.6 | 9.9 ± 0.6 |
| Real circ read, correct strand, pure boundary noise | 0.7 ± 0.4 | 0.0 ± 0.0 | 0.0 ± 0.0 | 43.8 ± 0.7 |
| Real circ read, but circRNA reported at a different/distant locus | 13.4 ± 0.7 | 20.6 ± 1.7 | 81.2 ± 7.6 | 46.3 ± 0.3 |
| Not from a circ read at all (linear read miscalled) | 0.0 ± 0.0 | 0.0 ± 0.0 | 0.0 ± 0.0 | 0.1 ± 0.1 |
| FP locus overlaps a genomic repeat | 19.8 ± 0.0 | 39.9 ± 0.7** | **55.2 ± 8.4** | **77.9 ± 0.2** |
| n FP loci (pooled, 3 runs) | 1,534 | 559 | 21 | 2,754 |
| n supporting reads (pooled, 3 runs) | 32,593 | 3,687 | 25 | 11,025 |


Many FP come from `antisense`-type circRNA, where the simulator deliberately sets the circRNA's strand opposite to its host gene. Simulator also mimics cDNA based RCRT-sequencing, meaning it simulates both strand sequences: one as an actual circRNA read and the other as cDNA read. The tools correctly find the genomic position of circRNA but select the host gene's (wrong) strand instead of flipping it to the antisense one (correct). Since strand agreement is required by all three matching strategies (exon-based, relaxed BSJ, strict BSJ), no boundary-tolerance setting can rescue these FP and we didn't find an easy way to override this gene-strand bias without introducing new FP elsewhere. These realistic strand misses could also happen on real data and we suggest to keep that in mind as you move onto wet-lab validation. 

**Practical implication**: most `discovery` false positives still point at the right gene/circRNA locus, just with an imprecise boundary or wrong strand. That makes them usable for gene- or locus-level expression signal analysis, but not for isoform-level analysis or primer/probe design, both of which need to get the exact strand/boundary right.

### What discovery still misses?

We also characterized every ground-truth circRNA against whether any tool recovered it at all (mean % ± SD across 3 human runs; missed + recovered = 100% of the GT catalog, ~8,700 circRNAs/run):

| | Missed (FN) | Recovered |
| --- | --- | --- |
| Share of GT catalog | 63.1 ± 0.1 | 36.9 ± 0.1 |
| Median mature length (bp) | 1,057 ± 12 | 291 ± 1 |
| Mean rolling-circle copy number | 7.6 ± 0.1 | 8.8 ± 0.0 |
| Median simulated read count | 34 ± 0 | 99 ± 1 |
| Canonical GT/AG splice site at BSJ | 63.0 ± 0.0 | 81.0 ± 0.0 |
| No gene overlap at all (intergenic) | 20.0 ± 0.0 | 9.0 ± 0.0 |
| Falls in a gene with ≥5 isoforms | 73.0 ± 0.0 | 84.0 ± 0.0 |
| Overlaps a genomic repeat | 51.0 ± 1.0 | 58.0 ± 1.0 |

Recall per biotype:

| Biotype | Missed | Recovered |
| --- | --- | --- |
| `eciRNA` (pure exonic) | 20.7 ± 0.2 | 79.3 ± 0.2 |
| `antisense` | 68.9 ± 0.1 | 31.1 ± 0.1 |
| `ciRNA` (purely intronic) | 74.5 ± 0.2 | 25.5 ± 0.2 |
| `EIciRNA` (exon+intron) | 76.4 ± 0.1 | 23.6 ± 0.1 |
| `intergenic` | 81.3 ± 0.2 | 18.7 ± 0.2 |

Missed circRNAs are generally longer (~3.6x), lower-expressed (~3x fewer reads), and use non-canonical splice sites more often. Rolling-circle copy number is also slightly lower, meaning weaker RCRT/RCA tandem-repeat evidence for the tools that rely on it. `eciRNA` seems to be the easier type to recover; every other biotype, especially without a clean match to an annotated exon, is more likely to be missed. Repeat overlap isn't a meaningful driver of circRNA discovery. 

### circRNA type classification

<img src="images/benchmark/type_classification_human_exon.png" width="100%"/>
<img src="images/benchmark/type_classification_human_relaxed_bsj.png" width="100%"/>
<img src="images/benchmark/type_classification_human_strict_bsj.png" width="100%"/>
<img src="images/benchmark/type_classification_mouse_exon.png" width="100%"/>
<img src="images/benchmark/type_classification_mouse_relaxed_bsj.png" width="100%"/>
<img src="images/benchmark/type_classification_mouse_strict_bsj.png" width="100%"/>

Same TP/FP-vs-GT breakdown as the GT coverage plot above, split by circRNA biotype (`eciRNA`, `EIciRNA`, `ciRNA`, `antisense`, `intergenic`, see [docs/output.md](output.md#circrna-type-classification)) instead of one bar per caller.

### Isoform detection

Long-read protocols are used mainly to resolve internal exon structure, not just BSJ position. Among ground-truth loci where the same BSJ has more than one real isoform (isoform family), we checked how often each caller recovers the exact isoform, versus calls the right BSJ with the wrong structure, versus misses the BSJ entirely.

<img src="images/benchmark/isoform_detection_reciprocal_overlap_99_human.png" width="100%"/>
<img src="images/benchmark/isoform_detection_reciprocal_overlap_99_mouse.png" width="100%"/>

<img src="images/benchmark/isoform_detection_reciprocal_overlap_90_human.png" width="100%"/>
<img src="images/benchmark/isoform_detection_reciprocal_overlap_90_mouse.png" width="100%"/>

Green+grey sum to 100% of the isoform-family ground truth; red shows the calls that found the right BSJ but didn't match internal structure criteria. 99% reciprocal exon overlap allows only trivial boundary noise; 90% allows a small sequence mismatch while still requiring the same structure, not a different one.

<img src="images/benchmark/isoform_family_completeness_reciprocal_overlap_99_human.png" width="100%"/>
<img src="images/benchmark/isoform_family_completeness_reciprocal_overlap_99_mouse.png" width="100%"/>

<img src="images/benchmark/isoform_family_completeness_reciprocal_overlap_90_human.png" width="100%"/>
<img src="images/benchmark/isoform_family_completeness_reciprocal_overlap_90_mouse.png" width="100%"/>

Per family plot instead of per isoform: did the caller get every isoform in the family right, some of them, or none. Getting one isoform right while missing or misclassifying other is the most common outcome; fully resolving isoform family is a rare event for all tools.

### Expression quantification and reproducibility

<img src="images/benchmark/expression_correlation_human.png" width="100%"/>
<img src="images/benchmark/expression_correlation_mouse.png" width="100%"/>

Read-count correlation against ground truth (right panel) and cross-method agreement (left panel), per tool/tier. Each cell shows three metrics:

- **Pearson r**: linear correlation between predicted and true counts. Sensitive to whether counts scale proportionally, not just in the same direction.
- **Spearman rho**: rank correlation. Only checks whether higher-expressed circRNAs get higher counts, regardless of exact scale. More resistant to outliers and non-linear (but still monotonic) relationships.
- **R²** (Pearson r squared): how much of the variance in true expression is captured by the predicted counts, on a 0-1 scale.

For differential expression, Spearman rho matters to rank circRNAs correctly relative to each other; Pearson r/R² matter to ensure the fold-changes DESeq2 estimates between conditions are quantitatively trustworthy, not just directionally right.

<img src="images/benchmark/sensitivity_by_expression_human.png" width="100%"/>
<img src="images/benchmark/sensitivity_by_expression_mouse.png" width="100%"/>

Detection sensitivity broken down by ground-truth expression tertile: low-expressed circRNAs are structurally harder to detect correctly than high-expressed ones, for every tool.

<img src="images/benchmark/reproducibility_pairs_human.png" width="100%"/>
<img src="images/benchmark/reproducibility_pairs_mouse.png" width="100%"/>

Pairwise read-count agreement across every tool/tier combination, ground truth included.

---

For how to run nanocirc and its command-line parameters, see the [nf_nanocirc_long_read pipeline](https://github.com/aerusakovich/nf_nanocirc_long_read).
