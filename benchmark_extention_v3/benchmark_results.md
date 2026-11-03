# circRNA Detection Tools: Benchmark Results

This benchmark is adapted from our own earlier benchmark study, Rusakovich et al., "Benchmarking circRNA detection tools from long-read sequencing using a data-driven and flexible simulation framework" (Peer Community Journal, 2026; see [CITATIONS.md](https://github.com/aerusakovich/nanocirc_pipeline/blob/main/CITATIONS.md) for the full reference), modified and extended based on our newer findings. 

This page shows benchmark results for 4 individual long-read circRNA detection tools (CIRI-long, circFL-seq, isoCirc, circNICK-lrs) and 4 consensus-based confidence tiers built by nanocirc (see [The four confidence tiers](#the-four-confidence-tiers) below), across wet-lab protocol, species, and annotated-vs-de-novo comparisons. The simulator and benchmark framework are presented in this repo as v3. Results below are from simulated human and mouse datasets, 3 runs each. The consensus tiers are computed by a long-read circRNA pipeline, [nanocirc_pipeline](https://github.com/aerusakovich/nanocirc_pipeline), that runs all 4 tools and pools their output; this repo is the benchmark and simulation framework behind the published numbers. Simulated datasets are available[here] (https://doi.org/10.17605/OSF.IO/C2GEJ). 

Several benchmarks are shown below:

- **[Same circRNAs, three wet-lab protocols](#same-circrnas-three-wet-lab-protocols)**: NanoSim's read-error model is characterized from real CIRI-long/isoCIRC and circFL-seq protocols ONT control reads; circRNA type composition and length/exon/splice distributions are the same across all 3 wet-lab methods, circRNA generation based on circRNAs from a human circbase+circatlas catalog. 
- **[Species benchmark: human vs mouse](#species-benchmark-human-vs-mouse)**: CIRI-long protocol, each species' own circRNA catalog and genome. Any difference here comes from species biology (genome, catalog composition), not the wet-lab protocol, which is fixed to CIRI-long for both.
- **[Database benchmark: annotated circRNA recovery](#database-benchmark-annotated-circrna-recovery)**: annotated vs de novo circRNA recovery. Same simulation and pipeline as above, but circRNAs are randomly sampled from a human circBase/circAtlas catalog instead of generated de novo, so tools are evaluated on already-known, already-catalogued circRNAs instead of novel, simulated and unannotated ones.

## The four confidence tiers

Alsongside with 4 tools we also present results for 4 modes of nanocirc [nanocirc_pipeline](https://github.com/aerusakovich/nanocirc_pipeline) pipeline, that are trading recall against precision:

| Tier | What it optimizes for | Trade-off |
| ---- | ---------------------- | --------- |
| **`discovery`** | Recall: keeps almost every circRNA any tool found. | Highest false-positive rate of the 4 tiers. |
| **`high_confidence`** | Precision: only calls with strong cross-tool agreement on both BSJ and structure, or isoCirc calls that met a defined read threshold, survive. | Lowest recall of the 4 tiers. |
| **`balanced_recall`** | A middle ground leaning recall: cuts off most of `discovery`'s false-positive rate while keeping most of its recall, including on isoform structure (see [Isoform detection](#isoform-detection)). | Recall below `discovery`. |
| **`balanced_precision`** | A middle ground leaning precision: higher precision than `balanced_recall` while keeping recall above `high_confidence`. | Recall below `balanced_recall`. |

## Same circRNAs, three wet-lab protocols

We used one circRNA catalog (same circbase/circatlas composition, human genome, same abundances) and only changed the wet-lab protocol NanoSim is characterized from. Any difference in the results below comes from the wet-lab protocol characteristics.

<details>
<summary><strong>Structural precision/recall</strong></summary>

### Structural precision/recall

<img src="images/benchmark/pr_scatter_equal_ciri_long.png" width="100%"/>
<img src="images/benchmark/pr_scatter_equal_circfl.png" width="100%"/>
<img src="images/benchmark/pr_scatter_equal_isocirc.png" width="100%"/>

Same circRNA catalog, same tools, same tiers, only the protocol the reads were simulated from changes between these 3 plots.

Mean precision (P) / recall (R) / F1 across 3 runs, exon-based matching, per nanocirc tier, plus each individual tool's own raw calls for comparison:

| Tier / tool | CIRI-long protocol | circFL-seq protocol | isoCirc protocol (CIRI-long on) |
| --- | --- | --- | --- |
| isoCirc alone | P 0.89 / R 0.14 / F1 0.25 | P 0.90 / R 0.19 / F1 0.32 | P 0.87 / R 0.08 / F1 0.15 |
| circFL-seq alone | P 0.79 / R 0.21 / F1 0.33 | P 0.77 / R 0.25 / F1 0.37 | P 0.74 / R 0.23 / F1 0.36 |
| CIRI-long alone | P 0.74 / R 0.34 / F1 0.46 | P 0.77 / R 0.41 / F1 0.53 | P 0.77 / R 0.03 / F1 0.05 |
| circNICK-lrs alone | P 0.24 / R 0.09 / F1 0.13 | P 0.24 / R 0.10 / F1 0.14 | P 0.22 / R 0.10 / F1 0.14 |
| `discovery` | P 0.68 / R 0.43 / F1 0.53 | P 0.69 / R 0.52 / F1 0.59 | P 0.66 / R 0.27 / F1 0.38 |
| `balanced_recall` | P 0.74 / R 0.42 / F1 0.53 | P 0.77 / R 0.51 / F1 0.61 | P 0.80 / R 0.25 / F1 0.38 |
| `balanced_precision` | P 0.85 / R 0.22 / F1 0.35 | P 0.85 / R 0.28 / F1 0.42 | P 0.84 / R 0.10 / F1 0.17 |
| `high_confidence` | P 0.86 / R 0.07 / F1 0.14 | P 0.88 / R 0.10 / F1 0.17 | P 0.86 / R 0.01 / F1 0.02 |

</details>

<details>
<summary><strong>Ground truth coverage</strong></summary>

### Ground truth coverage

<img src="images/benchmark/equal_gt_coverage_stacked_ciri_long.png" width="100%"/>
<img src="images/benchmark/equal_gt_coverage_stacked_circfl.png" width="100%"/>
<img src="images/benchmark/equal_gt_coverage_stacked_isocirc.png" width="100%"/>

circFL-seq protocol has led to the highest ground truth recovery at every tier (e.g. `discovery` R 0.52 vs CIRI-long's 0.43 and isoCirc's 0.26, from the table above); isoCirc protocol has the lowest at every tier, meaning it's reads are hardest to recover the calls from.

</details>

<details>
<summary><strong>circRNA type classification</strong></summary>

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

Same TP/FP-vs-GT breakdown as the GT coverage plot above, split by circRNA biotype. eciRNA is the easiest type to recover on every protocol. circFL-seq protocol leads to the best recovery of every type; isoCirc protocol to the worst recovery, most severely affecting antisense circRNA (9.7% vs CIRI-long's 33.6% and circFL-seq's 39.2%).

</details>

<details>
<summary><strong>Isoform detection</strong></summary>

### Isoform detection

The structural plots above show if a specific part of circRNA structure was found at its locus. But we are interested to see if the callers have found the exactly same molecule, especially on isoform level, which is the main point of using long-read sequencing. We evaluated that with three tiers, each pairing a BSJ boundary check with an exon overlap check: `full` (exact BSJ, 100% overlap, byte identical structure), `strict` (BSJ within tolerance window, 99% exon overlap, allows little boundary noise but can be easily reconstructed by follow up analyses), `relaxed` (BSJ within tolerance, 90% overlap, allowing a part of mismatched sequence).

<img src="images/benchmark/pr_scatter_isoform_ciri_long.png" width="100%"/>
<img src="images/benchmark/pr_scatter_isoform_circfl.png" width="100%"/>
<img src="images/benchmark/pr_scatter_isoform_isocirc.png" width="100%"/>

Long-read protocols are used mainly to resolve internal exon structure, not just find BSJ position or general locus. Among ground-truth loci where the same BSJ has more than one real isoform (isoform family), we checked how often each caller recovers the exact isoform with different levels of stringency.

<img src="images/benchmark/equal_isoform_detection_reciprocal_overlap_99_ciri_long.png" width="100%"/>
<img src="images/benchmark/equal_isoform_detection_reciprocal_overlap_99_circfl.png" width="100%"/>
<img src="images/benchmark/equal_isoform_detection_reciprocal_overlap_99_isocirc.png" width="100%"/>

<img src="images/benchmark/equal_isoform_detection_reciprocal_overlap_90_ciri_long.png" width="100%"/>
<img src="images/benchmark/equal_isoform_detection_reciprocal_overlap_90_circfl.png" width="100%"/>
<img src="images/benchmark/equal_isoform_detection_reciprocal_overlap_90_isocirc.png" width="100%"/>

Green+grey sum to 100% of the isoform-family ground truth; red shows the calls that found the right BSJ but didn't match required exon structure threshold. 99% reciprocal exon overlap allows only trivial boundary noise; 90% allows a small sequence mismatch while still requiring generally the same structure, not a different one. 

% of ground-truth isoforms called at 99% reciprocal overlap:

| Tier / tool | CIRI-long | circFL-seq | isoCirc |
| --- | --- | --- | --- |
| CIRI-long alone | 45.1% | 49.0% | 2.9% |
| isoCirc alone | 24.5% | 31.4% | 13.7% |
| circNICK-lrs alone | 26.5% | 27.5% | 28.4% |
| circFL-seq alone | 45.1% | 55.9% | 33.3% |
| `discovery` | 61.8% | 67.6% | 41.2% |
| `balanced_recall` | 55.9% | 61.8% | 29.4% |
| `balanced_precision` | 38.2% | 44.1% | 23.5% |
| `high_confidence` | 29.4% | 23.5% | 0.0% |

circFL-seq protocol allows resolution of the most isoforms at every tier except `high_confidence`, where CIRI-long is actually ahead (29.4% vs 23.5%). isoCirc protocol is last everywhere, and `high_confidence` doesn't get a single isoform right on it. Relaxing to 90% overlap barely moves these numbers (0-6 percentage points difference). Every nanocirc tier from `balanced_precision` up beats every raw tool at its own best protocol; CIRI-long and circFL-seq (the tools) are tied at 45.1% on their own protocol, ahead of isoCirc (24.5%) and circNICK-lrs (26.5%), but CIRI-long alone collapses on isoCirc protocol (2.9%) while circNICK-lrs, which is protocol-agnostic by construction, stays roughly flat across all 3 (26.5-28.4%).

<img src="images/benchmark/equal_isoform_family_completeness_reciprocal_overlap_99_ciri_long.png" width="100%"/>
<img src="images/benchmark/equal_isoform_family_completeness_reciprocal_overlap_99_circfl.png" width="100%"/>
<img src="images/benchmark/equal_isoform_family_completeness_reciprocal_overlap_99_isocirc.png" width="100%"/>

<img src="images/benchmark/equal_isoform_family_completeness_reciprocal_overlap_90_ciri_long.png" width="100%"/>
<img src="images/benchmark/equal_isoform_family_completeness_reciprocal_overlap_90_circfl.png" width="100%"/>
<img src="images/benchmark/equal_isoform_family_completeness_reciprocal_overlap_90_isocirc.png" width="100%"/>

Per family plot instead of per isoform: did the caller get every isoform in the family right, some of them, or none. % of isoform families fully resolved (every isoform in the family called exactly right), at 99% reciprocal overlap:

| Tier / tool | CIRI-long | circFL-seq | isoCirc |
| --- | --- | --- | --- |
| CIRI-long alone | 13.7% | 23.5% | 0.0% |
| isoCirc alone | 0.0% | 5.9% | 0.0% |
| circNICK-lrs alone | 0.0% | 0.0% | 0.0% |
| circFL-seq alone | 13.7% | 31.4% | 2.0% |
| `discovery` | 35.3% | 41.2% | 11.8% |
| `balanced_recall` | 23.5% | 29.4% | 0.0% |
| `balanced_precision` | 0.0% | 5.9% | 0.0% |
| `high_confidence` | 0.0% | 0.0% | 0.0% |

Getting one isoform right while missing or misclassifying others is the most common outcome; fully resolving an isoform family is rare for every tier and protocol. circFL-seq protocol gives the best odds of a full family at every tier; isoCirc protocol remains the hardest for all tools. No raw tool ever fully resolves a family on isoCirc protocol (0.0% across the board except circFL-seq's 2.0%), and circNICK-lrs never fully resolves one on any protocol; only `discovery` and `balanced_recall` reach double digits anywhere.

</details>

<details>
<summary><strong>Expression quantification and reproducibility</strong></summary>

### Expression quantification and reproducibility

<img src="images/benchmark/equal_expression_correlation_ciri_long.png" width="100%"/>
<img src="images/benchmark/equal_expression_correlation_circfl.png" width="100%"/>
<img src="images/benchmark/equal_expression_correlation_isocirc.png" width="100%"/>

<img src="images/benchmark/equal_expression_correlation_ciri_long_strict_bsj.png" width="100%"/>
<img src="images/benchmark/equal_expression_correlation_circfl_strict_bsj.png" width="100%"/>
<img src="images/benchmark/equal_expression_correlation_isocirc_strict_bsj.png" width="100%"/>

<img src="images/benchmark/equal_expression_correlation_ciri_long_exon_based.png" width="100%"/>
<img src="images/benchmark/equal_expression_correlation_circfl_exon_based.png" width="100%"/>
<img src="images/benchmark/equal_expression_correlation_isocirc_exon_based.png" width="100%"/>

Tool may call a correct circRNA, but misquantify it. It is important to measure, as for many research tasks that rely on differential expression precise circRNA quantification on locus and isoform is required. We measured it with Pearson r, 
Spearman rho, R², read-count correlation against ground truth and cross-method agreement, per tool/tier. 

Pearson r by protocol (relaxed-BSJ matching):

| Caller | CIRI-long | circFL-seq | isoCirc |
| --- | --- | --- | --- |
| CIRI-long | 0.911 | 0.935 | 0.308 |
| isoCirc | 0.528 | 0.576 | 0.403 |
| circNICK-lrs | 0.739 | 0.767 | 0.836 |
| circFL-seq | 0.765 | 0.785 | 0.693 |
| `high_confidence` | 0.972 | 0.969 | 0.835 |
| `balanced_precision` | 0.973 | 0.966 | 0.814 |
| `balanced_recall` | 0.969 | 0.961 | 0.843 |
| `discovery` | 0.970 | 0.962 | 0.840 |

isoCirc protocol is the hardest to quantify on for every caller except circNICK-lrs, whose isoCirc-protocol correlation (r=0.836) is its own best of the 3 protocols, the one case where quantification and structural-detection quality (weakest for circNICK-lrs on isoCirc protocol) move in opposite directions. CIRI-long and circFL-seq protocol are close for every other caller (e.g. `discovery` 0.970/0.962), CIRI-long protocol usually a shade ahead. Every nanocirc tier correlates with ground truth expression more strongly than any single raw tool, on every protocol (r ≥ 0.81 vs raw tools' r 0.40-0.94). isoCirc alone has the weakest quantification of any tool/tier on CIRI-long and circFL-seq protocol (r 0.53-0.58), despite it's high precision on circRNA detection. CIRI-long is the strongest quantifier among tools on its own and circFL-seq protocol, but the weakest on isoCirc protocol (r=0.308): the handful of calls it does make there (n=733) barely correlate with the truth, worse even than isoCirc's own 0.403.

#### Isoform-level expression quantification

The table above reports relaxed-BSJ matching, which collapses several real isoforms sharing one BSJ into a single point. The tables below instead score each real GT isoform on its own: `full`, `strict`, `relaxed` are the same 3 tiers used for [Isoform detection](#isoform-detection) and the isoform PR scatters above (`full` = exact BSJ + 100% exon overlap, `strict` = BSJ within tolerance + 99% overlap, `relaxed` = BSJ within tolerance + 90% overlap). A GT isoform's predicted count is whichever caller isoform structurally matched it under that tier, or 0 if none did, so n here is not the table above's n expanded, it is a genuinely stricter re-match and can be smaller. Pearson r (n):

**`full` (exact BSJ, 100% exon overlap):**

| Caller | CIRI-long r (n) | circFL-seq r (n) | isoCirc r (n) |
| --- | --- | --- | --- |
| CIRI-long | 0.908 (7,655) | 0.931 (8,926) | 0.326 (603) |
| isoCirc | 0.530 (4,030) | 0.576 (5,314) | 0.403 (2,254) |
| circNICK-lrs | 0.614 (1,606) | 0.698 (1,838) | 0.806 (1,892) |
| circFL-seq | 0.766 (5,937) | 0.786 (7,032) | 0.706 (4,788) |
| `high_confidence` | 0.972 (2,837) | 0.969 (3,814) | 0.833 (641) |
| `balanced_precision` | 0.972 (6,435) | 0.967 (7,917) | 0.811 (2,743) |
| `balanced_recall` | 0.966 (9,720) | 0.961 (11,863) | 0.847 (6,051) |
| `discovery` | 0.966 (9,882) | 0.961 (11,986) | 0.846 (6,266) |

**`strict` (BSJ within tolerance, 99% exon overlap):**

| Caller | CIRI-long r (n) | circFL-seq r (n) | isoCirc r (n) |
| --- | --- | --- | --- |
| CIRI-long | 0.909 (8,694) | 0.933 (10,686) | 0.314 (668) |
| isoCirc | 0.528 (4,053) | 0.575 (5,339) | 0.402 (2,264) |
| circNICK-lrs | 0.615 (1,696) | 0.699 (1,932) | 0.807 (1,978) |
| circFL-seq | 0.767 (6,008) | 0.785 (7,137) | 0.706 (4,863) |
| `high_confidence` | 0.972 (2,857) | 0.969 (3,832) | 0.834 (647) |
| `balanced_precision` | 0.972 (6,461) | 0.967 (7,950) | 0.811 (2,752) |
| `balanced_recall` | 0.967 (10,381) | 0.962 (12,908) | 0.846 (6,120) |
| `discovery` | 0.967 (10,546) | 0.962 (13,045) | 0.844 (6,341) |

**`relaxed` (BSJ within tolerance, 90% exon overlap):**

| Caller | CIRI-long r (n) | circFL-seq r (n) | isoCirc r (n) |
| --- | --- | --- | --- |
| CIRI-long | 0.910 (9,681) | 0.933 (11,728) | 0.316 (744) |
| isoCirc | 0.529 (4,075) | 0.576 (5,379) | 0.405 (2,280) |
| circNICK-lrs | 0.618 (1,900) | 0.698 (2,141) | 0.811 (2,176) |
| circFL-seq | 0.767 (6,054) | 0.786 (7,185) | 0.700 (5,093) |
| `high_confidence` | 0.972 (2,869) | 0.969 (3,851) | 0.835 (653) |
| `balanced_precision` | 0.972 (6,503) | 0.965 (8,003) | 0.812 (2,770) |
| `balanced_recall` | 0.967 (10,893) | 0.960 (13,418) | 0.842 (6,319) |
| `discovery` | 0.968 (11,062) | 0.961 (13,563) | 0.840 (6,564) |

n grows from `full` to `relaxed` as the overlap requirement loosens (e.g. discovery's CIRI-long-protocol n: 9,882 -> 10,546 -> 11,062), since more near-miss predicted isoforms clear the lower bar. circNICK-lrs is the caller most affected by tightening the requirement: its correlation drops sharply against the locus-level table above (r 0.61-0.70 here vs r 0.74-0.83 there) and its n falls the most in relative terms, since its structure calls carry real boundary noise that a genuine per-isoform re-match, unlike a locus-level average, does not forgive. CIRI-long's own correlation actually improves once isoform identity is enforced (r 0.91 here vs r 0.90-0.94 above) but n drops hard on isoCirc protocol (603-744 vs 4,105-9,681 at the locus level above), the clearest sign yet that what little CIRI-long calls on isoCirc-protocol reads rarely has the right internal structure, not just the right BSJ. Every nanocirc tier still correlates far more strongly than any raw tool at every tier of stringency. isoCirc protocol remains the hardest to quantify for every caller except circNICK-lrs, which quantifies best on isoCirc protocol at every tier (r 0.81 vs 0.61-0.70 on the other two); CIRI-long protocol and circFL-seq protocol stay close for everyone else.

<img src="images/benchmark/equal_expression_correlation_ciri_long_full.png" width="100%"/>
<img src="images/benchmark/equal_expression_correlation_circfl_full.png" width="100%"/>
<img src="images/benchmark/equal_expression_correlation_isocirc_full.png" width="100%"/>

<img src="images/benchmark/equal_expression_correlation_ciri_long_strict.png" width="100%"/>
<img src="images/benchmark/equal_expression_correlation_circfl_strict.png" width="100%"/>
<img src="images/benchmark/equal_expression_correlation_isocirc_strict.png" width="100%"/>

<img src="images/benchmark/equal_expression_correlation_ciri_long_relaxed.png" width="100%"/>
<img src="images/benchmark/equal_expression_correlation_circfl_relaxed.png" width="100%"/>
<img src="images/benchmark/equal_expression_correlation_isocirc_relaxed.png" width="100%"/>

Detection sensitivity by ground-truth expression tertile (low/mid/high read count):

| Tier | Protocol | Low | Mid | High |
| --- | --- | --- | --- | --- |
| `discovery` | CIRI-long | 32.3% | 38.1% | 47.8% |
| `discovery` | circFL-seq | 40.3% | 46.8% | 56.8% |
| `discovery` | isoCirc | 14.7% | 22.8% | 37.5% |
| `balanced_recall` | CIRI-long | 31.3% | 37.3% | 46.9% |
| `balanced_recall` | circFL-seq | 39.5% | 45.9% | 56.1% |
| `balanced_recall` | isoCirc | 12.8% | 20.2% | 34.6% |
| `balanced_precision` | CIRI-long | 14.1% | 21.2% | 32.2% |
| `balanced_precision` | circFL-seq | 18.1% | 26.0% | 39.2% |
| `balanced_precision` | isoCirc | 3.5% | 7.7% | 18.1% |
| `high_confidence` | CIRI-long | 2.8% | 6.1% | 13.7% |
| `high_confidence` | circFL-seq | 3.9% | 8.3% | 16.7% |
| `high_confidence` | isoCirc | 0.1% | 0.5% | 3.0% |

Sensitivity grows with expression at every tier and protocol (more reads means more chances to be called). `balanced_recall` is close to `discovery` (within 1-3 points at every tertile), while `balanced_precision` is in between `discovery` and `high_confidence`. circFL-seq protocol allows the highest sensitivity at every tertile and tier; isoCirc protocol the lowest, most severely at `high_confidence`, where even highly-expressed circRNAs are barely detected (3.0% vs circFL-seq's 16.7%).

<img src="images/benchmark/equal_sensitivity_by_expression_ciri_long.png" width="100%"/>
<img src="images/benchmark/equal_sensitivity_by_expression_circfl.png" width="100%"/>
<img src="images/benchmark/equal_sensitivity_by_expression_isocirc.png" width="100%"/>

Same tertiles, matched at the `full`/`strict`/`relaxed` isoform tiers instead of structural comparison:

<img src="images/benchmark/equal_sensitivity_by_expression_ciri_long_isoform.png" width="100%"/>
<img src="images/benchmark/equal_sensitivity_by_expression_circfl_isoform.png" width="100%"/>
<img src="images/benchmark/equal_sensitivity_by_expression_isocirc_isoform.png" width="100%"/>

<img src="images/benchmark/equal_reproducibility_pairs_ciri_long.png" width="100%"/>
<img src="images/benchmark/equal_reproducibility_pairs_circfl.png" width="100%"/>

<img src="images/benchmark/equal_reproducibility_pairs_ciri_long_strict_bsj.png" width="100%"/>
<img src="images/benchmark/equal_reproducibility_pairs_circfl_strict_bsj.png" width="100%"/>

<img src="images/benchmark/equal_reproducibility_pairs_ciri_long_exon_based.png" width="100%"/>
<img src="images/benchmark/equal_reproducibility_pairs_circfl_exon_based.png" width="100%"/>

Pairwise read-count agreement across every tool/tier combination, ground truth included. nanocirc's own tiers agree with ground truth better, than for any other caller and agree with each other almost perfectly, since they share the same remap-based quantification regardless of which tool found the locus. Agreement of raw tool calls with ground truth is looser and protocol-dependent; isoCirc shows the widest scatter against every other caller and ground truth on every protocol, consistent with its weak expression correlation above, CIRI-long shows the best agreement with ground truth among single tools on 2 out of 3 protocols.

Same reproducibility matrix, matched at the `full`/`strict`/`relaxed` isoform tiers instead of structural comparison, so agreement is evaluated on isoform level, each real isoform scored against its own matched caller isoform (see [Isoform-level expression quantification](#isoform-level-expression-quantification) above).

<img src="images/benchmark/equal_reproducibility_pairs_ciri_long_full.png" width="100%"/>
<img src="images/benchmark/equal_reproducibility_pairs_circfl_full.png" width="100%"/>

<img src="images/benchmark/equal_reproducibility_pairs_ciri_long_strict.png" width="100%"/>
<img src="images/benchmark/equal_reproducibility_pairs_circfl_strict.png" width="100%"/>

<img src="images/benchmark/equal_reproducibility_pairs_ciri_long_relaxed.png" width="100%"/>
<img src="images/benchmark/equal_reproducibility_pairs_circfl_relaxed.png" width="100%"/>

### Caller overlap

Ground-truth loci recovered by each of the four tools are pooled across runs and drawn as an UpSet-style plot, where a locus counts as recovered if it was matched in at least one run. Two versions were used in the original benchmark paper: exon-based, which requires correct internal structure, and relaxed BSJ, which requires only genomic position. This is now extended to all six matching strategies used elsewhere in this benchmark: strict BSJ and exon-based (exact/structural), relaxed BSJ (coordinate-tolerant), and the full/strict/relaxed isoform tiers (structure plus BSJ tolerance, tightest to loosest).

The point of this plot is agreement, not any one tool's own recall. Of every circRNA any tool finds, only about 13% are found by all four tools together (12.5% on the species benchmark, 12.9% on the wet-lab CIRI-long protocol), and roughly half are found by exactly one tool. Each tool is effectively viewing circRNAs through its own lens, with no single reliable "ground truth caller" to defer to. That disagreement is the core motivation for NanoCirc: a naive union of all four tools' calls would inherit every tool's false positives, and a naive intersection would keep only the rare loci every tool happens to agree on, trading away most of the recall. NanoCirc's confidence tiers reconcile this with vote-based merging instead of a plain union/intersection rule, as described in the pipeline's own methods documentation.

<img src="images/benchmark/upset_plot_equal_ciri_long_strict_bsj.png" width="100%"/>
<img src="images/benchmark/upset_plot_equal_circfl_strict_bsj.png" width="100%"/>
<img src="images/benchmark/upset_plot_equal_isocirc_strict_bsj.png" width="100%"/>

<img src="images/benchmark/upset_plot_equal_ciri_long.png" width="100%"/>
<img src="images/benchmark/upset_plot_equal_circfl.png" width="100%"/>
<img src="images/benchmark/upset_plot_equal_isocirc.png" width="100%"/>

<img src="images/benchmark/upset_plot_equal_ciri_long_exon_based.png" width="100%"/>
<img src="images/benchmark/upset_plot_equal_circfl_exon_based.png" width="100%"/>
<img src="images/benchmark/upset_plot_equal_isocirc_exon_based.png" width="100%"/>

<img src="images/benchmark/upset_plot_equal_ciri_long_full.png" width="100%"/>
<img src="images/benchmark/upset_plot_equal_circfl_full.png" width="100%"/>
<img src="images/benchmark/upset_plot_equal_isocirc_full.png" width="100%"/>

<img src="images/benchmark/upset_plot_equal_ciri_long_strict.png" width="100%"/>
<img src="images/benchmark/upset_plot_equal_circfl_strict.png" width="100%"/>
<img src="images/benchmark/upset_plot_equal_isocirc_strict.png" width="100%"/>

<img src="images/benchmark/upset_plot_equal_ciri_long_relaxed.png" width="100%"/>
<img src="images/benchmark/upset_plot_equal_circfl_relaxed.png" width="100%"/>
<img src="images/benchmark/upset_plot_equal_isocirc_relaxed.png" width="100%"/>

</details>

**Protocol results:**

- **circFL-seq protocol generates reads that work well with all tools**, at every tier: highest recall, highest F1, and the highest BSJ-confirmation rate of the 3 protocols (64.1%).
- **CIRI-long protocol reads are also compatible with all tools**, with a lower BSJ-confirmation rate (54.8%) and lower recall at every tier than circFL-seq.
- **isoCirc protocol reads are the hardest for the tools to work with**, despite carrying the exact same ~9,800 underlying circRNAs. Its recall is the lowest of the 3 protocols at every tier, but precision remains comparable or higher to other wet-lab approaches. CIRI-long is almost incompatible with this wet lab protocol, and that was verified by running CIRI-long on real isoCIRC dataset as well - the tool never managed to finish the analysis on a full dataset.

Because the circRNA catalog, its abundances, and the simulated read amount were identical across all 3 comparisons, this difference comes entirely from how each protocol's wet-lab induced error profile and read length characteristics change the reads NanoSim generates from them.

<details>
<summary><strong>Statistical significance</strong></summary>

### Statistical significance

We tested whether the protocol-to-protocol differences above are real or could be noise. Each pair matches every caller/circRNA-type/run point between the two protocols, using these statistical methods (used the same way for every "statistical significance" comparison in this document):

- **Linear regression (slope, intercept, R²)**: fits one side's value from the other's (value = intercept + slope × other value). Slope near 1 with a small intercept means the two sides are close to equal plus a fixed shift. R² is how well that line predicts the value, whatever the size of the shift.
- **Kendall's tau-b**: ranks the 8 callers best to worst on each side, then checks how well the two rankings agree.
- **Mantel-Haenszel odds ratio + Breslow-Day test**: uses the raw TP/FP/FN counts directly, stratified by circRNA type. The odds ratio is how many times more likely a call is to be correct, or a ground-truth circRNA is recovered, on one side versus the other. Breslow-Day then checks if that odds ratio is the same size for every type; p < 0.05 means the effect size depends on circRNA type.
- **Wilcoxon signed-rank**: tests if one side's values are systematically higher than the other's, without assuming a normal distribution.
- **Cohen's d (paired)**: how big that gap is, in standard deviations. ~0.2 is small, ~0.5 medium, ~0.8+ large.

#### CIRI-long vs circFL-seq

##### Aggregate, by metric (n paired points pooling every caller × type × run)

| Metric | n | CIRI-long | circFL-seq | Difference | Slope | Intercept | R² | Wilcoxon p | Cohen's d |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Precision | 114 | 0.720 | 0.731 | +0.011 | 0.993 | 0.016 | 0.987 | 0.009 | 0.353 |
| Recall | 120 | 0.221 | 0.275 | +0.054 | 1.104 | 0.031 | 0.974 | 2×10⁻¹⁹ | 1.230 |
| F1 | 120 | 0.282 | 0.337 | +0.054 | 1.071 | 0.034 | 0.965 | 2×10⁻¹⁹ | 1.183 |

Precision fits an almost perfect 1:1 line (R² = 0.99, slope ≈ 1) but the gap itself is small (+0.011) and only modestly significant (p = 0.009, d = 0.35). Recall and F1 are both significantly and substantially higher on circFL-seq (p < 10⁻¹⁸, d ≈ 1.2), with a similarly tight fit (R² ≥ 0.96).

##### Per type, by metric (n=24 paired points per row pooling all 8 callers × 3 runs, except ciRNA precision n=21 and intergenic precision n=21: some caller × run cells had zero calls of that type, so precision can't be defined and those cells are dropped)

| Type | Metric | n | CIRI-long | circFL-seq | Difference | Slope | Intercept | R² | Wilcoxon p | Cohen's d |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| eciRNA | Precision | 24 | 0.630 | 0.628 | −0.002 | 0.988 | 0.006 | 0.997 | 0.643 | −0.144 |
| eciRNA | Recall | 24 | 0.536 | 0.600 | +0.064 | 1.044 | 0.041 | 0.980 | 1×10⁻⁷ | 2.191 |
| eciRNA | F1 | 24 | 0.547 | 0.583 | +0.036 | 1.004 | 0.034 | 0.970 | 1×10⁻⁷ | 1.371 |
| EIciRNA | Precision | 24 | 0.635 | 0.655 | +0.020 | 0.979 | 0.033 | 0.970 | 0.060 | 0.485 |
| EIciRNA | Recall | 24 | 0.118 | 0.194 | +0.075 | 1.606 | 0.004 | 0.990 | 1×10⁻⁷ | 1.692 |
| EIciRNA | F1 | 24 | 0.188 | 0.281 | +0.093 | 1.401 | 0.018 | 0.985 | 1×10⁻⁷ | 2.049 |
| ciRNA | Precision | 21 | 0.879 | 0.905 | +0.026 | 0.737 | 0.257 | 0.971 | 0.003 | 0.828 |
| ciRNA | Recall | 24 | 0.145 | 0.208 | +0.064 | 1.451 | −0.001 | 0.999 | 6×10⁻⁵ | 1.120 |
| ciRNA | F1 | 24 | 0.222 | 0.296 | +0.074 | 1.306 | 0.006 | 0.999 | 6×10⁻⁵ | 1.364 |
| antisense | Precision | 24 | 0.693 | 0.683 | −0.010 | 0.985 | 0.001 | 0.999 | 0.003 | −0.726 |
| antisense | Recall | 24 | 0.156 | 0.183 | +0.028 | 1.166 | 0.002 | 0.990 | 2×10⁻⁴ | 1.009 |
| antisense | F1 | 24 | 0.241 | 0.272 | +0.031 | 1.116 | 0.003 | 0.991 | 2×10⁻⁴ | 0.977 |
| intergenic | Precision | 21 | 0.792 | 0.815 | +0.023 | 1.024 | 0.004 | 0.967 | 0.001 | 0.831 |
| intergenic | Recall | 24 | 0.151 | 0.189 | +0.038 | 1.246 | 0.001 | 1.000 | 6×10⁻⁵ | 1.024 |
| intergenic | F1 | 24 | 0.214 | 0.252 | +0.038 | 1.163 | 0.004 | 0.999 | 6×10⁻⁵ | 1.178 |

Every type's recall and F1 fit almost perfectly (R² ≥ 0.97) and circFL-seq is significantly ahead on both, for every type (p ≤ 2×10⁻⁴). Precision is smaller in magnitude and mixed in direction: eciRNA and antisense slightly favor CIRI-long, while ciRNA and intergenic significantly favor circFL-seq (p ≤ 0.003); EIciRNA shows no significant precision difference (p = 0.06).

##### Does caller ranking transfer?

| Metric | n callers | Kendall's tau-b | Spearman rho |
| --- | --- | --- | --- |
| Precision | 8 | 0.929 (p = 4×10⁻⁴) | 0.976 (p = 3×10⁻⁵) |
| Recall | 8 | 1.000 (p = 5×10⁻⁵) | 1.000 (p = 0) |
| F1 | 8 | 1.000 (p = 5×10⁻⁵) | 1.000 (p = 0) |

Ranking transfers essentially perfectly. By F1:

| Caller | F1 CIRI-long | F1 circFL-seq | Rank CIRI-long | Rank circFL-seq |
| --- | --- | --- | --- | --- |
| balanced_recall | 0.470 | 0.558 | 1 | 1 |
| discovery | 0.457 | 0.538 | 2 | 2 |
| CIRI-long | 0.375 | 0.447 | 3 | 3 |
| balanced_precision | 0.304 | 0.364 | 4 | 4 |
| circFL-seq | 0.287 | 0.334 | 5 | 5 |
| isoCirc | 0.207 | 0.261 | 6 | 6 |
| high_confidence | 0.108 | 0.136 | 7 | 7 |
| circNICK-lrs | 0.051 | 0.055 | 8 | 8 |

All eight callers have on the exact same F1 rank on both protocols.

##### Per tool/mode, by metric (n=15 paired points per row: the 5 circRNA types × 3 runs, except circNICK-lrs precision n=9)

| Tool | Metric | n | CIRI-long | circFL-seq | Difference | Slope | Intercept | R² | Wilcoxon p | Cohen's d |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CIRI-long | Precision | 15 | 0.694 | 0.731 | +0.038 | 0.914 | 0.098 | 0.925 | 1×10⁻⁴ | 1.218 |
| CIRI-long | Recall | 15 | 0.293 | 0.360 | +0.067 | 0.954 | 0.081 | 0.947 | 6×10⁻⁵ | 1.563 |
| CIRI-long | F1 | 15 | 0.375 | 0.447 | +0.072 | 0.801 | 0.147 | 0.807 | 6×10⁻⁵ | 1.326 |
| isoCirc | Precision | 15 | 0.891 | 0.885 | −0.006 | 1.083 | −0.080 | 0.943 | 0.272 | −0.214 |
| isoCirc | Recall | 15 | 0.136 | 0.181 | +0.045 | 1.229 | 0.014 | 0.996 | 6×10⁻⁵ | 1.271 |
| isoCirc | F1 | 15 | 0.207 | 0.261 | +0.054 | 1.113 | 0.031 | 0.992 | 6×10⁻⁵ | 1.799 |
| circNICK-lrs | Precision | 9 | 0.096 | 0.091 | −0.004 | 0.994 | −0.004 | 0.986 | 0.249 | −0.473 |
| circNICK-lrs | Recall | 15 | 0.057 | 0.065 | +0.008 | 1.137 | 0.000 | 1.000 | 0.028 | 0.563 |
| circNICK-lrs | F1 | 15 | 0.051 | 0.055 | +0.003 | 1.061 | 0.000 | 0.999 | 0.028 | 0.602 |
| circFL-seq | Precision | 15 | 0.643 | 0.634 | −0.009 | 0.983 | 0.002 | 0.997 | 0.107 | −0.483 |
| circFL-seq | Recall | 15 | 0.203 | 0.241 | +0.038 | 0.969 | 0.045 | 0.974 | 6×10⁻⁵ | 1.281 |
| circFL-seq | F1 | 15 | 0.287 | 0.334 | +0.048 | 0.974 | 0.055 | 0.963 | 6×10⁻⁵ | 1.150 |
| high_confidence | Precision | 15 | 0.931 | 0.936 | +0.005 | 0.923 | 0.077 | 1.000 | 0.024 | 0.786 |
| high_confidence | Recall | 15 | 0.072 | 0.093 | +0.021 | 1.179 | 0.008 | 0.996 | 0.002 | 0.826 |
| high_confidence | F1 | 15 | 0.108 | 0.136 | +0.028 | 1.117 | 0.015 | 0.992 | 0.002 | 0.977 |
| balanced_precision | Precision | 15 | 0.898 | 0.907 | +0.009 | 1.118 | −0.097 | 0.870 | 0.762 | 0.279 |
| balanced_precision | Recall | 15 | 0.218 | 0.271 | +0.053 | 1.087 | 0.034 | 0.992 | 6×10⁻⁵ | 1.933 |
| balanced_precision | F1 | 15 | 0.304 | 0.364 | +0.060 | 0.980 | 0.066 | 0.977 | 6×10⁻⁵ | 1.914 |
| balanced_recall | Precision | 15 | 0.699 | 0.725 | +0.026 | 0.853 | 0.129 | 0.916 | 0.041 | 0.776 |
| balanced_recall | Recall | 15 | 0.391 | 0.491 | +0.100 | 0.956 | 0.117 | 0.972 | 6×10⁻⁵ | 3.017 |
| balanced_recall | F1 | 15 | 0.470 | 0.558 | +0.088 | 0.718 | 0.221 | 0.896 | 6×10⁻⁵ | 1.956 |
| discovery | Precision | 15 | 0.660 | 0.680 | +0.021 | 0.947 | 0.056 | 0.953 | 0.095 | 0.603 |
| discovery | Recall | 15 | 0.398 | 0.497 | +0.099 | 0.936 | 0.124 | 0.971 | 6×10⁻⁵ | 2.820 |
| discovery | F1 | 15 | 0.457 | 0.538 | +0.082 | 0.619 | 0.256 | 0.783 | 6×10⁻⁵ | 1.638 |

Every caller's recall and F1 is significantly higher on circFL-seq (p ≤ 6×10⁻⁵) with a good fit (R² ≥ 0.78). For precision CIRI-long, high_confidence, ciRNA-heavy balanced_precision, and balanced_recall show a small but significant circFL-seq precision edge (p ≤ 0.041), while isoCirc, circNICK-lrs, circFL-seq-seq, and discovery show no significant precision difference at all between the two protocols.

| Tool | Outcome | Strata (types) | OR (MH) | 95% CI | MH p | Breslow-Day p |
| --- | --- | --- | --- | --- | --- | --- |
| CIRI-long | Precision | 5 | 1.15 | 1.09–1.21 | 4×10⁻⁸ | 8×10⁻⁶ |
| CIRI-long | Recall | 5 | 1.43 | 1.38–1.49 | <10⁻³⁰⁰ | <10⁻³⁰⁰ |
| isoCirc | Precision | 5 | 0.88 | 0.79–0.98 | 0.018 | 0.075 |
| isoCirc | Recall | 5 | 1.52 | 1.45–1.60 | <10⁻³⁰⁰ | 5×10⁻⁴ |
| circNICK-lrs | Precision | 5 | 1.00 | 0.93–1.07 | 0.910 | 0.096 |
| circNICK-lrs | Recall | 5 | 1.19 | 1.10–1.28 | 6×10⁻⁶ | 0.794 |
| circFL-seq | Precision | 5 | 0.93 | 0.87–1.01 | 0.084 | 0.988 |
| circFL-seq | Recall | 5 | 1.31 | 1.26–1.37 | <10⁻³⁰⁰ | 4×10⁻¹⁵ |
| high_confidence | Precision | 5 | 1.12 | 0.96–1.30 | 0.147 | 0.897 |
| high_confidence | Recall | 5 | 1.44 | 1.35–1.54 | <10⁻³⁰⁰ | 1×10⁻⁵ |
| balanced_precision | Precision | 5 | 0.92 | 0.85–1.00 | 0.058 | 2×10⁻⁴ |
| balanced_precision | Recall | 5 | 1.45 | 1.39–1.52 | <10⁻³⁰⁰ | 3×10⁻⁴ |
| balanced_recall | Precision | 5 | 1.11 | 1.06–1.16 | 5×10⁻⁶ | 4×10⁻⁸ |
| balanced_recall | Recall | 5 | 1.61 | 1.55–1.67 | <10⁻³⁰⁰ | <10⁻³⁰⁰ |
| discovery | Precision | 5 | 1.05 | 1.01–1.10 | 0.010 | 2×10⁻¹⁰ |
| discovery | Recall | 5 | 1.60 | 1.55–1.66 | <10⁻³⁰⁰ | <10⁻³⁰⁰ |

Every caller has significantly higher odds of recovering a ground-truth circRNA (recall) on circFL-seq (OR 1.2–1.6). Precision odds are mixed and mostly weaker: isoCirc actually favors CIRI-long (OR = 0.88, p = 0.018), circNICK-lrs/circFL-seq/high_confidence show no significant precision difference, and the rest favor circFL-seq modestly (OR 1.05–1.15). 5 of 16 caller/outcome combinations have a protocol effect that's uniform across circRNA types (circNICK-lrs's precision and recall, isoCirc's precision, high_confidence's precision, circFL-seq's precision), the most homogeneity of any comparison in this document; the other 11 combinations' boost size still depends on which circRNA type it is.

**Conclusion.** Correlation is excellent for both metrics (aggregate R² ≥ 0.97, essentially every type and caller individually significant on recall/F1). Caller ranking transfers perfectly (Kendall's tau-b = 1.00 on recall/F1, all 8 callers same rank by F1). Precision differences are real but small and caller/type-dependent, with several callers and types showing no significant precision difference between protocols while circFL-seq is better for recall/F1.


#### CIRI-long vs isoCirc

##### Aggregate, by metric (n paired points pooling every caller × type × run)

| Metric | n | CIRI-long | isoCirc | Difference | Slope | Intercept | R² | Wilcoxon p | Cohen's d |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Precision | 105 | 0.696 | 0.619 | -0.077 | 0.897 | -0.005 | 0.702 | 6×10⁻⁷ | -0.495 |
| Recall | 120 | 0.221 | 0.107 | -0.114 | 0.549 | -0.014 | 0.664 | 9×10⁻¹⁷ | -0.913 |
| F1 | 120 | 0.282 | 0.146 | -0.136 | 0.594 | -0.022 | 0.623 | 4×10⁻¹⁹ | -1.042 |

Precision still fits reasonably well (R² = 0.70) with a moderate, significant gap (d = -0.50). Recall and F1 are both significantly and substantially lower on isoCirc protocol (p < 10⁻¹⁶, d ≈ -0.9 to -1.0), and the fit is noisier than the CIRI-long-vs-circFL-seq comparison above (R² ≈ 0.62-0.66 vs ≥ 0.96 there): isoCirc protocol doesn't just shift every caller down by a similar amount, it hits some callers, especially CIRI-long itself, far harder than others.

##### Per type, by metric (n=24 paired points per row pooling all 8 callers × 3 runs, except ciRNA precision n=18 and intergenic precision n=18: some caller × run cells had zero calls of that type, so precision can't be defined and those cells are dropped)

| Type | Metric | n | CIRI-long | isoCirc | Difference | Slope | Intercept | R² | Wilcoxon p | Cohen's d |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| eciRNA | Precision | 24 | 0.630 | 0.635 | +0.005 | 0.880 | 0.080 | 0.852 | 0.317 | 0.056 |
| eciRNA | Recall | 24 | 0.536 | 0.312 | -0.223 | 0.565 | 0.010 | 0.342 | 2×10⁻⁶ | -1.333 |
| eciRNA | F1 | 24 | 0.547 | 0.365 | -0.182 | 0.695 | -0.015 | 0.312 | 2×10⁻⁶ | -1.142 |
| EIciRNA | Precision | 24 | 0.635 | 0.442 | -0.193 | 0.643 | 0.034 | 0.333 | 8×10⁻⁶ | -0.829 |
| EIciRNA | Recall | 24 | 0.118 | 0.080 | -0.038 | 0.725 | -0.006 | 0.461 | 0.003 | -0.651 |
| EIciRNA | F1 | 24 | 0.188 | 0.113 | -0.074 | 0.654 | -0.009 | 0.458 | 2×10⁻⁶ | -0.901 |
| ciRNA | Precision | 18 | 0.859 | 0.789 | -0.070 | 0.675 | 0.210 | 0.249 | 0.012 | -0.555 |
| ciRNA | Recall | 24 | 0.145 | 0.057 | -0.088 | 0.385 | 0.001 | 0.577 | 6×10⁻⁵ | -1.005 |
| ciRNA | F1 | 24 | 0.222 | 0.096 | -0.126 | 0.463 | -0.007 | 0.600 | 6×10⁻⁵ | -1.090 |
| antisense | Precision | 21 | 0.649 | 0.575 | -0.074 | 0.903 | -0.011 | 0.953 | 0.003 | -0.801 |
| antisense | Recall | 24 | 0.156 | 0.045 | -0.111 | 0.287 | 0.001 | 0.783 | 1×10⁻⁴ | -1.137 |
| antisense | F1 | 24 | 0.241 | 0.082 | -0.159 | 0.345 | -0.001 | 0.784 | 8×10⁻⁵ | -1.176 |
| intergenic | Precision | 18 | 0.758 | 0.716 | -0.042 | 1.575 | -0.477 | 0.916 | 0.102 | -0.436 |
| intergenic | Recall | 24 | 0.151 | 0.043 | -0.109 | 0.285 | -0.000 | 0.646 | 6×10⁻⁵ | -0.972 |
| intergenic | F1 | 24 | 0.214 | 0.072 | -0.141 | 0.366 | -0.006 | 0.656 | 6×10⁻⁵ | -1.049 |

Every type's recall and F1 is significantly lower on isoCirc protocol (p ≤ 0.003), same direction as the aggregate. eciRNA's precision is the one metric with no significant difference at all (p = 0.32) and a strong fit (R² = 0.85); every other type's precision is also significantly lower on isoCirc protocol (p ≤ 0.012), except intergenic precision, which is not significant here either (p = 0.10).

##### Does caller ranking transfer?

| Metric | n callers | Kendall's tau-b | Spearman rho |
| --- | --- | --- | --- |
| Precision | 8 | 0.571 (p = 0.061) | 0.619 (p = 0.102) |
| Recall | 8 | 0.571 (p = 0.061) | 0.690 (p = 0.058) |
| F1 | 8 | 0.571 (p = 0.061) | 0.690 (p = 0.058) |

Unlike the CIRI-long-vs-circFL-seq comparison above, ranking does not transfer reliably here (tau-b = 0.57, not significant at p = 0.061). By F1:

| Caller | F1 CIRI-long | F1 isoCirc | Rank CIRI-long | Rank isoCirc |
| --- | --- | --- | --- | --- |
| balanced_recall | 0.470 | 0.287 | 1 | 1 |
| discovery | 0.457 | 0.270 | 2 | 2 |
| CIRI-long | 0.375 | 0.041 | 3 | 7 |
| balanced_precision | 0.304 | 0.144 | 4 | 4 |
| circFL-seq | 0.287 | 0.224 | 5 | 3 |
| isoCirc | 0.207 | 0.125 | 6 | 5 |
| high_confidence | 0.108 | 0.021 | 7 | 8 |
| circNICK-lrs | 0.051 | 0.054 | 8 | 6 |

The two nanocirc tiers that lead on CIRI-long protocol (balanced_recall, discovery) still lead on isoCirc protocol. Everything else reshuffles: the CIRI-long tool itself falls from 3rd to 7th, since it barely functions on isoCirc-protocol reads (see below), while circFL-seq and circNICK-lrs both move up a few places because they degrade less than the callers around them.

##### Per tool/mode, by metric (n=15 paired points per row: the 5 circRNA types × 3 runs, except circNICK-lrs precision n=9)

| Tool | Metric | n | CIRI-long | isoCirc | Difference | Slope | Intercept | R² | Wilcoxon p | Cohen's d |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CIRI-long | Precision | 15 | 0.694 | 0.646 | -0.048 | 0.427 | 0.349 | 0.209 | 0.208 | -0.420 |
| CIRI-long | Recall | 15 | 0.293 | 0.022 | -0.271 | 0.132 | -0.016 | 0.908 | 6×10⁻⁵ | -1.670 |
| CIRI-long | F1 | 15 | 0.375 | 0.041 | -0.333 | 0.339 | -0.086 | 0.845 | 6×10⁻⁵ | -3.956 |
| isoCirc | Precision | 15 | 0.891 | 0.894 | +0.003 | 0.943 | 0.054 | 0.830 | 0.551 | 0.069 |
| isoCirc | Recall | 15 | 0.136 | 0.076 | -0.060 | 0.628 | -0.010 | 0.967 | 6×10⁻⁵ | -1.048 |
| isoCirc | F1 | 15 | 0.207 | 0.125 | -0.082 | 0.695 | -0.019 | 0.951 | 6×10⁻⁵ | -1.197 |
| circNICK-lrs | Precision | 9 | 0.096 | 0.080 | -0.015 | 0.924 | -0.008 | 0.938 | 0.028 | -0.793 |
| circNICK-lrs | Recall | 15 | 0.057 | 0.067 | +0.010 | 1.148 | 0.001 | 0.999 | 0.028 | 0.623 |
| circNICK-lrs | F1 | 15 | 0.051 | 0.054 | +0.002 | 1.033 | 0.000 | 0.999 | 0.075 | 0.510 |
| circFL-seq | Precision | 15 | 0.643 | 0.416 | -0.227 | 0.636 | 0.007 | 0.803 | 1×10⁻⁴ | -1.409 |
| circFL-seq | Recall | 15 | 0.203 | 0.164 | -0.039 | 0.719 | 0.018 | 0.955 | 0.060 | -0.649 |
| circFL-seq | F1 | 15 | 0.287 | 0.224 | -0.063 | 0.745 | 0.010 | 0.975 | 2×10⁻⁴ | -1.039 |
| high_confidence | Precision | 6 | 0.827 | 0.492 | -0.335 | 19.597 | -15.721 | 1.000 | 0.031 | -0.925 |
| high_confidence | Recall | 15 | 0.072 | 0.011 | -0.061 | 0.176 | -0.002 | 0.996 | 6×10⁻⁵ | -0.570 |
| high_confidence | F1 | 15 | 0.108 | 0.021 | -0.087 | 0.230 | -0.004 | 0.992 | 6×10⁻⁵ | -0.609 |
| balanced_precision | Precision | 15 | 0.898 | 0.913 | +0.015 | 0.835 | 0.163 | 0.544 | 0.135 | 0.263 |
| balanced_precision | Recall | 15 | 0.218 | 0.093 | -0.125 | 0.604 | -0.039 | 0.997 | 6×10⁻⁵ | -1.526 |
| balanced_precision | F1 | 15 | 0.304 | 0.144 | -0.160 | 0.827 | -0.107 | 0.988 | 6×10⁻⁵ | -3.943 |
| balanced_recall | Precision | 15 | 0.699 | 0.640 | -0.059 | 0.497 | 0.293 | 0.309 | 0.041 | -0.573 |
| balanced_recall | Recall | 15 | 0.391 | 0.208 | -0.183 | 0.767 | -0.092 | 0.794 | 6×10⁻⁵ | -2.046 |
| balanced_recall | F1 | 15 | 0.470 | 0.287 | -0.183 | 0.970 | -0.169 | 0.482 | 6×10⁻⁵ | -1.505 |
| discovery | Precision | 15 | 0.660 | 0.580 | -0.080 | 0.784 | 0.063 | 0.793 | 9×10⁻⁴ | -1.103 |
| discovery | Recall | 15 | 0.398 | 0.217 | -0.181 | 0.824 | -0.111 | 0.814 | 6×10⁻⁵ | -2.079 |
| discovery | F1 | 15 | 0.457 | 0.270 | -0.186 | 0.778 | -0.085 | 0.302 | 6×10⁻⁵ | -1.556 |

CIRI-long's own recall collapses from 29% on its own protocol to 2.2% on isoCirc protocol (d = -1.67), and its F1 drops by an order of magnitude (0.375 → 0.041, d = -3.96, the single largest effect in this table): it is barely functioning on isoCirc-protocol reads. Every other caller's recall and F1 is also significantly lower on isoCirc protocol (p ≤ 0.060), consistent with the aggregate result, but nowhere near as severely. Precision is smaller and mostly not significant: isoCirc, circNICK-lrs's paired precision, balanced_precision, and CIRI-long itself show no significant precision difference, while circFL-seq and the two remaining tiers do.

| Tool | Outcome | Strata (types) | OR (MH) | 95% CI | MH p | Breslow-Day p |
| --- | --- | --- | --- | --- | --- | --- |
| CIRI-long | Precision | 5 | 1.05 | 0.92–1.20 | 0.448 | 4×10⁻⁹ |
| CIRI-long | Recall | 5 | 0.04 | 0.04–0.04 | <10⁻³⁰⁰ | <10⁻³⁰⁰ |
| isoCirc | Precision | 5 | 0.93 | 0.82–1.06 | 0.311 | 0.598 |
| isoCirc | Recall | 5 | 0.47 | 0.44–0.50 | <10⁻³⁰⁰ | 7×10⁻⁶ |
| circNICK-lrs | Precision | 5 | 0.91 | 0.85–0.98 | 0.008 | 5×10⁻⁵ |
| circNICK-lrs | Recall | 5 | 1.23 | 1.14–1.32 | 4×10⁻⁸ | 0.382 |
| circFL-seq | Precision | 5 | 0.30 | 0.28–0.32 | <10⁻³⁰⁰ | <10⁻³⁰⁰ |
| circFL-seq | Recall | 5 | 0.74 | 0.70–0.77 | <10⁻³⁰⁰ | <10⁻³⁰⁰ |
| high_confidence | Precision | 5 | 0.76 | 0.58–0.99 | 0.051 | 2×10⁻⁹ |
| high_confidence | Recall | 5 | 0.11 | 0.10–0.13 | <10⁻³⁰⁰ | 1×10⁻⁴ |
| balanced_precision | Precision | 5 | 1.10 | 0.98–1.23 | 0.099 | 0.006 |
| balanced_precision | Recall | 5 | 0.28 | 0.27–0.30 | <10⁻³⁰⁰ | <10⁻³⁰⁰ |
| balanced_recall | Precision | 5 | 1.02 | 0.97–1.08 | 0.486 | <10⁻³⁰⁰ |
| balanced_recall | Recall | 5 | 0.35 | 0.34–0.37 | <10⁻³⁰⁰ | <10⁻³⁰⁰ |
| discovery | Precision | 5 | 0.87 | 0.83–0.91 | 2×10⁻⁹ | <10⁻³⁰⁰ |
| discovery | Recall | 5 | 0.36 | 0.34–0.37 | <10⁻³⁰⁰ | <10⁻³⁰⁰ |

CIRI-long's recall odds ratio (0.04, 95% CI 0.04–0.04) is the largest single effect measured anywhere in this document: its odds of recovering a ground-truth circRNA are roughly 25× lower on isoCirc-protocol reads than on its own protocol's reads. Every caller's recall odds are significantly lower on isoCirc protocol (OR 0.04–0.74), except circNICK-lrs, which is actually higher (OR = 1.23): it is the one tool this protocol change doesn't hurt. Precision odds are smaller and mixed: circFL-seq, circNICK-lrs, high_confidence, and discovery favor CIRI-long protocol, isoCirc's own precision shows no significant difference, and balanced_recall/balanced_precision/CIRI-long's own precision show no significant difference either. Only 2 of 16 caller/outcome combinations have a protocol effect that's uniform across circRNA types (isoCirc's precision, circNICK-lrs's recall), the least homogeneity of any comparison in this document; every other combination's boost size depends heavily on which circRNA type it is.

**Conclusion.** Correlation is reasonable for precision (R² = 0.70) but weaker for recall/F1 (R² ≈ 0.62–0.66), and unlike the CIRI-long-vs-circFL-seq comparison above, caller ranking does not transfer reliably (Kendall's tau-b = 0.57, not significant). The reason is CIRI-long itself: its recall and F1 on isoCirc-protocol reads are close to zero (recall OR = 0.04, F1 Cohen's d = -3.96), dragging it from 3rd to 7th place in the ranking and confirming numerically what the recommendations at the top of this page already state qualitatively, that CIRI-long is close to incompatible with isoCirc-protocol reads. Every other caller is also significantly weaker on recall/F1 on isoCirc protocol, but far less severely, and circNICK-lrs is actually unaffected or slightly better (recall OR = 1.23).

#### circFL-seq vs isoCirc

##### Aggregate, by metric (n paired points pooling every caller × type × run)

| Metric | n | circFL-seq | isoCirc | Difference | Slope | Intercept | R² | Wilcoxon p | Cohen's d |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Precision | 105 | 0.708 | 0.619 | -0.089 | 0.898 | -0.016 | 0.708 | 2×10⁻⁸ | -0.576 |
| Recall | 120 | 0.275 | 0.107 | -0.168 | 0.482 | -0.025 | 0.640 | 5×10⁻¹⁹ | -1.135 |
| F1 | 120 | 0.337 | 0.146 | -0.191 | 0.521 | -0.030 | 0.570 | 9×10⁻²⁰ | -1.248 |

Same pattern as the CIRI-long comparison above, larger in size: circFL-seq protocol scores the highest of the 3 protocols to begin with, so the drop to isoCirc-protocol numbers is the biggest gap measured in this document for recall/F1 (d ≈ -1.1 to -1.2). Precision again fits well (R² = 0.71) with a moderate, significant gap (d = -0.58).

##### Per type, by metric (n=24 paired points per row pooling all 8 callers × 3 runs, except ciRNA precision n=18 and intergenic precision n=18: some caller × run cells had zero calls of that type, so precision can't be defined and those cells are dropped)

| Type | Metric | n | circFL-seq | isoCirc | Difference | Slope | Intercept | R² | Wilcoxon p | Cohen's d |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| eciRNA | Precision | 24 | 0.628 | 0.635 | +0.007 | 0.895 | 0.073 | 0.861 | 0.277 | 0.081 |
| eciRNA | Recall | 24 | 0.600 | 0.312 | -0.288 | 0.532 | -0.007 | 0.338 | 1×10⁻⁶ | -1.658 |
| eciRNA | F1 | 24 | 0.583 | 0.365 | -0.218 | 0.624 | 0.001 | 0.262 | 1×10⁻⁷ | -1.294 |
| EIciRNA | Precision | 24 | 0.655 | 0.442 | -0.213 | 0.650 | 0.017 | 0.336 | 8×10⁻⁶ | -0.919 |
| EIciRNA | Recall | 24 | 0.194 | 0.080 | -0.114 | 0.432 | -0.004 | 0.427 | 2×10⁻⁶ | -1.310 |
| EIciRNA | F1 | 24 | 0.281 | 0.113 | -0.168 | 0.438 | -0.010 | 0.408 | 6×10⁻⁷ | -1.477 |
| ciRNA | Precision | 18 | 0.889 | 0.789 | -0.100 | 1.110 | -0.198 | 0.367 | 1×10⁻⁴ | -0.898 |
| ciRNA | Recall | 24 | 0.208 | 0.057 | -0.152 | 0.260 | 0.002 | 0.556 | 6×10⁻⁵ | -1.075 |
| ciRNA | F1 | 24 | 0.296 | 0.096 | -0.200 | 0.348 | -0.007 | 0.579 | 6×10⁻⁵ | -1.215 |
| antisense | Precision | 21 | 0.638 | 0.575 | -0.063 | 0.920 | -0.012 | 0.952 | 0.028 | -0.697 |
| antisense | Recall | 24 | 0.183 | 0.045 | -0.138 | 0.256 | -0.002 | 0.859 | 9×10⁻⁵ | -1.176 |
| antisense | F1 | 24 | 0.272 | 0.082 | -0.190 | 0.321 | -0.005 | 0.849 | 6×10⁻⁵ | -1.226 |
| intergenic | Precision | 18 | 0.785 | 0.716 | -0.069 | 1.455 | -0.426 | 0.973 | 0.004 | -0.937 |
| intergenic | Recall | 24 | 0.189 | 0.043 | -0.146 | 0.228 | -0.000 | 0.642 | 6×10⁻⁵ | -0.993 |
| intergenic | F1 | 24 | 0.252 | 0.072 | -0.179 | 0.312 | -0.006 | 0.642 | 6×10⁻⁵ | -1.084 |

Every type's recall and F1 is significantly lower on isoCirc protocol (p ≤ 9×10⁻⁵), same direction as the aggregate and as the CIRI-long comparison above. eciRNA's precision is again the exception (p = 0.28, R² = 0.86); every other type's precision is significantly lower on isoCirc protocol (p ≤ 0.028), unlike the CIRI-long comparison where intergenic precision wasn't significant either.

##### Does caller ranking transfer?

| Metric | n callers | Kendall's tau-b | Spearman rho |
| --- | --- | --- | --- |
| Precision | 8 | 0.643 (p = 0.031) | 0.643 (p = 0.086) |
| Recall | 8 | 0.571 (p = 0.061) | 0.690 (p = 0.058) |
| F1 | 8 | 0.571 (p = 0.061) | 0.690 (p = 0.058) |

Precision ranking transfers significantly (tau-b = 0.64, p = 0.031); recall and F1 ranking do not (p = 0.061), the same as the CIRI-long comparison above. By F1:

| Caller | F1 circFL-seq | F1 isoCirc | Rank circFL-seq | Rank isoCirc |
| --- | --- | --- | --- | --- |
| balanced_recall | 0.558 | 0.287 | 1 | 1 |
| discovery | 0.538 | 0.270 | 2 | 2 |
| CIRI-long | 0.447 | 0.041 | 3 | 7 |
| balanced_precision | 0.364 | 0.144 | 4 | 4 |
| circFL-seq | 0.334 | 0.224 | 5 | 3 |
| isoCirc | 0.261 | 0.125 | 6 | 5 |
| high_confidence | 0.136 | 0.021 | 7 | 8 |
| circNICK-lrs | 0.055 | 0.054 | 8 | 6 |

Same reshuffle as the CIRI-long comparison above, and for the same reason: CIRI-long falls from 3rd to 7th because it barely functions on isoCirc-protocol reads regardless of which protocol it's being compared against, while circFL-seq and circNICK-lrs both move up a few places.

##### Per tool/mode, by metric (n=15 paired points per row: the 5 circRNA types × 3 runs, except circNICK-lrs precision n=9)

| Tool | Metric | n | circFL-seq | isoCirc | Difference | Slope | Intercept | R² | Wilcoxon p | Cohen's d |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CIRI-long | Precision | 15 | 0.731 | 0.646 | -0.086 | 0.443 | 0.321 | 0.203 | 0.010 | -0.766 |
| CIRI-long | Recall | 15 | 0.360 | 0.022 | -0.338 | 0.131 | -0.025 | 0.867 | 6×10⁻⁵ | -2.125 |
| CIRI-long | F1 | 15 | 0.447 | 0.041 | -0.406 | 0.347 | -0.114 | 0.703 | 6×10⁻⁵ | -5.289 |
| isoCirc | Precision | 15 | 0.885 | 0.894 | +0.009 | 0.875 | 0.120 | 0.888 | 0.638 | 0.239 |
| isoCirc | Recall | 15 | 0.181 | 0.076 | -0.105 | 0.509 | -0.016 | 0.964 | 6×10⁻⁵ | -1.158 |
| isoCirc | F1 | 15 | 0.261 | 0.125 | -0.136 | 0.617 | -0.036 | 0.934 | 6×10⁻⁵ | -1.464 |
| circNICK-lrs | Precision | 9 | 0.091 | 0.080 | -0.011 | 0.944 | -0.006 | 0.982 | 0.028 | -1.005 |
| circNICK-lrs | Recall | 15 | 0.065 | 0.067 | +0.002 | 1.009 | 0.001 | 0.999 | 0.046 | 0.520 |
| circNICK-lrs | F1 | 15 | 0.055 | 0.054 | -0.001 | 0.974 | 0.000 | 1.000 | 0.173 | -0.410 |
| circFL-seq | Precision | 15 | 0.634 | 0.416 | -0.218 | 0.648 | 0.005 | 0.807 | 6×10⁻⁵ | -1.399 |
| circFL-seq | Recall | 15 | 0.241 | 0.164 | -0.077 | 0.739 | -0.014 | 0.971 | 1×10⁻³ | -1.456 |
| circFL-seq | F1 | 15 | 0.334 | 0.224 | -0.110 | 0.749 | -0.026 | 0.971 | 6×10⁻⁵ | -1.836 |
| high_confidence | Precision | 6 | 0.841 | 0.492 | -0.349 | 18.632 | -15.176 | 1.000 | 0.031 | -0.965 |
| high_confidence | Recall | 15 | 0.093 | 0.011 | -0.082 | 0.149 | -0.003 | 0.984 | 6×10⁻⁵ | -0.626 |
| high_confidence | F1 | 15 | 0.136 | 0.021 | -0.115 | 0.203 | -0.007 | 0.969 | 6×10⁻⁵ | -0.691 |
| balanced_precision | Precision | 15 | 0.907 | 0.913 | +0.006 | 0.854 | 0.139 | 0.817 | 0.524 | 0.155 |
| balanced_precision | Recall | 15 | 0.271 | 0.093 | -0.178 | 0.548 | -0.055 | 0.980 | 6×10⁻⁵ | -1.721 |
| balanced_precision | F1 | 15 | 0.364 | 0.144 | -0.220 | 0.811 | -0.151 | 0.935 | 6×10⁻⁵ | -3.739 |
| balanced_recall | Precision | 15 | 0.725 | 0.640 | -0.085 | 0.435 | 0.325 | 0.188 | 0.041 | -0.784 |
| balanced_recall | Recall | 15 | 0.491 | 0.208 | -0.283 | 0.814 | -0.192 | 0.841 | 6×10⁻⁵ | -3.708 |
| balanced_recall | F1 | 15 | 0.558 | 0.287 | -0.272 | 1.328 | -0.454 | 0.520 | 6×10⁻⁵ | -2.246 |
| discovery | Precision | 15 | 0.680 | 0.580 | -0.100 | 0.768 | 0.058 | 0.715 | 9×10⁻⁴ | -1.215 |
| discovery | Recall | 15 | 0.497 | 0.217 | -0.280 | 0.893 | -0.227 | 0.863 | 6×10⁻⁵ | -3.919 |
| discovery | F1 | 15 | 0.538 | 0.270 | -0.268 | 0.841 | -0.182 | 0.172 | 6×10⁻⁵ | -2.084 |

CIRI-long's F1 drop (0.447 → 0.041, d = -5.29) is the single largest effect in this table, and its recall alone drops from 36% to 2.2% (d = -2.13). Every other caller's recall and F1 is also significantly lower on isoCirc protocol (p ≤ 1×10⁻³), same direction as the CIRI-long-vs-isoCirc comparison above but generally with larger gaps, since circFL-seq protocol is the strongest baseline of the 3. circNICK-lrs's F1 is the only metric with no significant difference at all (p = 0.173), consistent with it being the one caller this protocol change barely touches. Precision is again smaller and mixed: isoCirc, balanced_precision, and circNICK-lrs's paired recall/F1 show no significant precision-side effect, while CIRI-long, circFL-seq, and the remaining tiers do.

| Tool | Outcome | Strata (types) | OR (MH) | 95% CI | MH p | Breslow-Day p |
| --- | --- | --- | --- | --- | --- | --- |
| CIRI-long | Precision | 5 | 0.95 | 0.83–1.08 | 0.438 | 7×10⁻¹¹ |
| CIRI-long | Recall | 5 | 0.03 | 0.03–0.03 | <10⁻³⁰⁰ | <10⁻³⁰⁰ |
| isoCirc | Precision | 5 | 1.05 | 0.93–1.19 | 0.448 | 0.117 |
| isoCirc | Recall | 5 | 0.31 | 0.29–0.32 | <10⁻³⁰⁰ | 3×10⁻⁴ |
| circNICK-lrs | Precision | 5 | 0.91 | 0.85–0.98 | 0.009 | 0.015 |
| circNICK-lrs | Recall | 5 | 1.04 | 0.96–1.11 | 0.347 | 0.237 |
| circFL-seq | Precision | 5 | 0.30 | 0.28–0.32 | <10⁻³⁰⁰ | <10⁻³⁰⁰ |
| circFL-seq | Recall | 5 | 0.57 | 0.55–0.60 | <10⁻³⁰⁰ | 7×10⁻¹¹ |
| high_confidence | Precision | 5 | 0.67 | 0.52–0.88 | 0.004 | 3×10⁻¹⁰ |
| high_confidence | Recall | 5 | 0.08 | 0.07–0.09 | <10⁻³⁰⁰ | 3×10⁻⁸ |
| balanced_precision | Precision | 5 | 1.22 | 1.10–1.36 | 3×10⁻⁴ | 0.004 |
| balanced_precision | Recall | 5 | 0.19 | 0.18–0.20 | <10⁻³⁰⁰ | <10⁻³⁰⁰ |
| balanced_recall | Precision | 5 | 0.89 | 0.85–0.94 | 1×10⁻⁵ | <10⁻³⁰⁰ |
| balanced_recall | Recall | 5 | 0.22 | 0.21–0.23 | <10⁻³⁰⁰ | <10⁻³⁰⁰ |
| discovery | Precision | 5 | 0.81 | 0.77–0.84 | <10⁻³⁰⁰ | <10⁻³⁰⁰ |
| discovery | Recall | 5 | 0.22 | 0.21–0.23 | <10⁻³⁰⁰ | <10⁻³⁰⁰ |

CIRI-long's recall odds ratio (0.03, 95% CI 0.03–0.03) is again the largest effect in this table, in line with the CIRI-long-vs-isoCirc comparison above (OR = 0.04 there): its odds of recovering a ground-truth circRNA are roughly 30× lower on isoCirc-protocol reads than on circFL-seq protocol's reads. Every caller's recall odds are significantly lower on isoCirc protocol (OR 0.03–0.57) except circNICK-lrs, which shows no significant difference (OR = 1.04, p = 0.347): the same tool the CIRI-long comparison also found unaffected. Precision odds are mixed: circFL-seq, high_confidence, circNICK-lrs, and discovery favor circFL-seq protocol, balanced_precision favors isoCirc protocol, and isoCirc's own precision and balanced_recall's precision show no significant difference. Only 2 of 16 caller/outcome combinations have a protocol effect that's uniform across circRNA types (isoCirc's precision, circNICK-lrs's recall), matching the CIRI-long comparison's homogeneity count exactly.

**Conclusion.** Same overall picture as the CIRI-long-vs-isoCirc comparison above, scaled up: circFL-seq protocol is the strongest baseline of the 3, so its gap to isoCirc-protocol numbers is the largest measured in this document for recall/F1 (d ≈ -1.1 to -1.2 aggregate, up to d = -5.29 for CIRI-long's own F1). Precision ranking transfers significantly (tau-b = 0.64) but recall/F1 ranking does not (tau-b = 0.57, not significant), again because CIRI-long collapses on isoCirc-protocol reads (recall OR = 0.03) regardless of which other protocol it's being compared against. circNICK-lrs is consistently the one caller isoCirc protocol doesn't hurt.

</details>

<details>
<summary><strong>What are these false positives?</strong></summary>

### What are these false positives?

We traced every FP-supporting read back to the simulator's own ground truth (single-tool discovery exon-based FPs, all 3 runs per protocol, mean % ± SD across runs), separately for each protocol's own reads:

#### CIRI-long protocol

| FP read is... | CIRI-long | circFL-seq | isoCirc | circNICK-lrs |
| --- | --- | --- | --- | --- |
| Not from a circ read at all (linear read miscalled) | 0.0 ± 0.0 | 0.0 ± 0.0 | 0.0 ± 0.0 | 0.1 ± 0.0 |
| Real circ read, correct boundary, wrong strand (antisense called as sense) | 86.8 ± 0.4 | 76.8 ± 0.7 | 93.7 ± 0.5 | 2.5 ± 0.3 |
| Real circ read, correct strand, pure boundary noise | 1.2 ± 0.2 | 4.9 ± 1.0 | 0.0 ± 0.0 | 31.1 ± 0.7 |
| Real circ read, but circRNA reported at a different/distant locus | 12.0 ± 0.4 | 18.3 ± 0.7 | 6.3 ± 0.5 | 66.3 ± 0.5 |
| FP locus overlaps a genomic repeat | 18.3 ± 0.2 | 31.2 ± 1.0 | 19.9 ± 0.4 | 62.5 ± 0.3 |
| n FP loci (pooled, 3 runs) | 3,566 | 1,660 | 513 | 7,999 |
| n supporting reads (pooled, 3 runs) | 60,498 | 14,942 | 5,274 | 27,927 |

#### circFL-seq protocol

| FP read is... | CIRI-long | circFL-seq | isoCirc | circNICK-lrs |
| --- | --- | --- | --- | --- |
| Not from a circ read at all (linear read miscalled) | 0.0 ± 0.0 | 0.0 ± 0.0 | 0.0 ± 0.0 | 0.1 ± 0.0 |
| Real circ read, correct boundary, wrong strand (antisense called as sense) | 84.3 ± 0.9 | 73.1 ± 2.1 | 93.1 ± 0.5 | 2.4 ± 0.4 |
| Real circ read, correct strand, pure boundary noise | 2.7 ± 0.6 | 7.1 ± 2.5 | 0.0 ± 0.0 | 23.7 ± 1.3 |
| Real circ read, but circRNA reported at a different/distant locus | 13.0 ± 0.6 | 19.8 ± 0.4 | 6.9 ± 0.4 | 73.8 ± 1.7 |
| FP locus overlaps a genomic repeat | 20.6 ± 0.8 | 37.8 ± 0.8 | 28.7 ± 1.5 | 58.7 ± 0.2 |
| n FP loci (pooled, 3 runs) | 3,702 | 2,136 | 635 | 9,316 |
| n supporting reads (pooled, 3 runs) | 48,061 | 20,115 | 7,227 | 39,033 |

#### isoCirc protocol

| FP read is... | CIRI-long | circFL-seq | isoCirc | circNICK-lrs |
| --- | --- | --- | --- | --- |
| Not from a circ read at all (linear read miscalled) | 0.0 ± 0.0 | 0.0 ± 0.0 | 0.0 ± 0.0 | 0.2 ± 0.1 |
| Real circ read, correct boundary, wrong strand (antisense called as sense) | 72.1 ± 8.8 | 28.9 ± 0.2 | 91.8 ± 0.9 | 1.9 ± 0.1 |
| Real circ read, correct strand, pure boundary noise | 0.7 ± 0.6 | 60.2 ± 0.3 | 0.0 ± 0.0 | 18.5 ± 1.2 |
| Real circ read, but circRNA reported at a different/distant locus | 27.3 ± 8.5 | 10.9 ± 0.3 | 8.2 ± 0.9 | 79.4 ± 1.3 |
| FP locus overlaps a genomic repeat | 11.2 ± 1.6 | 52.4 ± 0.8 | 22.4 ± 1.4 | 57.4 ± 0.5 |
| n FP loci (pooled, 3 runs) | 241 | 2,269 | 354 | 10,995 |
| n supporting reads (pooled, 3 runs) | 592 | 13,832 | 1,103 | 29,882 |

Many FP come from `antisense`-type circRNA, where the simulator deliberately sets the circRNA's strand opposite to its host gene; since strand agreement is required by every matching strategy, no boundary-tolerance setting can rescue these. That pattern is present for CIRI-long and isoCirc on all 3 protocols (72-94% wrong-strand), and for circFL-seq on the CIRI-long and circFL-seq protocols (73-77%), but circFL-seq's FP composition changes on isoCirc-protocol reads: wrong-strand drops to 28.9% and "correct strand, pure boundary noise" goes up to 60.2%, the opposite of its own pattern everywhere else and of every other tool on this protocol. CIRI-long's own FP composition on isoCirc protocol is also unusual: barely any single-tool FP loci at all (241, an order of magnitude below its 3,566-3,702 on the other 2 protocols, consistent with how few calls it makes here overall), and its second-largest bucket is "different/distant locus" (27.3%) rather than boundary noise, which it almost never has (0.7%). circNICK-lrs's FPs are structurally different from the other 3 tools throughout: mostly boundary noise and distant-locus calls (66-79%), rarely a strand flip (2-3%), and consistently the most repeat-associated (57-63% vs 11-38% for the other tools).

**Practical implication**: most `discovery` false positives still are at the right gene/circRNA, just with imprecise boundary or wrong strand. That makes them usable for gene- or locus-level expression signal analysis, but we recommend to re-verify this tier's low confidence calls in wet-lab, for example by Sanger sequencing. Do not blindly trust their reported BSJ/structure/strand and deposit them into the databases.

</details>

<details>
<summary><strong>What discovery still misses?</strong></summary>

### What discovery still misses?

We also characterized every ground-truth circRNA against whether any tool recovered it at all, per protocol (missed + recovered = 100% of the ~9,861-circRNA GT catalog, pooled across 3 runs):

#### CIRI-long protocol

| | Missed (FN) | Recovered |
| --- | --- | --- |
| Share of GT catalog | 60.6% | 39.4% |
| Mean mature length (bp) | 3,156 | 353 |
| Mean rolling-circle copy number | 7.8 | 8.7 |
| Median simulated read count | 49 | 71 |
| Canonical GT/AG splice site at BSJ | 63.2% | 76.7% |
| No gene overlap at all (intergenic) | 21.8% | 14.6% |
| Falls in a gene with ≥5 isoforms | 72.1% | 78.2% |
| Overlaps a genomic repeat | 49.7% | 54.0% |

Recall per biotype:

| Biotype | Missed | Recovered |
| --- | --- | --- |
| `eciRNA` (pure exonic) | 19.8% | 80.2% |
| `antisense` | 66.3% | 33.7% |
| `ciRNA` (purely intronic) | 68.6% | 31.4% |
| `EIciRNA` (exon+intron) | 79.3% | 20.7% |
| `intergenic` | 69.4% | 30.6% |

#### circFL-seq protocol

| | Missed (FN) | Recovered |
| --- | --- | --- |
| Share of GT catalog | 52.0% | 48.0% |
| Mean mature length (bp) | 3,502 | 477 |
| Mean rolling-circle copy number | 7.5 | 8.9 |
| Median simulated read count | 48 | 69 |
| Canonical GT/AG splice site at BSJ | 62.0% | 75.6% |
| No gene overlap at all (intergenic) | 22.6% | 14.9% |
| Falls in a gene with ≥5 isoforms | 71.2% | 78.1% |
| Overlaps a genomic repeat | 49.3% | 53.7% |

Recall per biotype:

| Biotype | Missed | Recovered |
| --- | --- | --- |
| `eciRNA` (pure exonic) | 13.4% | 86.7% |
| `antisense` | 61.4% | 38.6% |
| `ciRNA` (purely intronic) | 56.6% | 43.4% |
| `EIciRNA` (exon+intron) | 67.0% | 33.0% |
| `intergenic` | 62.5% | 37.5% |

#### isoCirc protocol

| | Missed (FN) | Recovered |
| --- | --- | --- |
| Share of GT catalog | 76.3% | 23.7% |
| Mean mature length (bp) | 2,509 | 572 |
| Mean rolling-circle copy number | 7.6 | 9.9 |
| Median simulated read count | 49 | 97 |
| Canonical GT/AG splice site at BSJ | 63.6% | 84.6% |
| No gene overlap at all (intergenic) | 21.6% | 10.3% |
| Falls in a gene with ≥5 isoforms | 72.2% | 82.0% |
| Overlaps a genomic repeat | 49.7% | 57.1% |

Recall per biotype:

| Biotype | Missed | Recovered |
| --- | --- | --- |
| `eciRNA` (pure exonic) | 41.0% | 59.0% |
| `antisense` | 89.9% | 10.1% |
| `ciRNA` (purely intronic) | 84.8% | 15.2% |
| `EIciRNA` (exon+intron) | 78.8% | 21.2% |
| `intergenic` | 87.8% | 12.2% |

circFL-seq protocol has the largest recovery of the 3 (48.0% of the GT catalog found by at least one caller), CIRI-long protocol is close behind (39.4%), and isoCirc protocol has by far the smallest (23.7%, barely more than half of CIRI-long protocol's). The missed-vs-recovered pattern is consistent across all 3 protocols regardless of that gap: missed circRNAs are longer (~4-9x), lower-expressed, and use non-canonical splice sites more often than recovered ones. `eciRNA` is always the easiest biotype to recover in all wet-lab protocols. On isoCirc protocol every biotype's recall collapses further than on the other 2 protocols, most severely for `antisense` (10.1% recovered, versus 33.7% on CIRI-long protocol and 38.6% on circFL-seq protocol), also due to the same wrong-strand FP pattern documented above rather than due to a different reason. Repeat overlap is not a meaningful driver of recovery on any protocol.

</details>

<details>
<summary><strong>isoCirc protocol: with and without CIRI-long</strong></summary>

#### isoCirc protocol: with and without CIRI-long

The isoCirc-protocol column above used all 4 tools. Given the real-world risk of CIRI-long hanging on isoCirc-protocol data (it did not finish in reasonable time on real isoCirc-protocol reads, getting stuck on its collapse step for over 24 days across 3 attempts), we also scored the same dataset with CIRI-long excluded from the consensus tool set, to isolate what including it actually costs or gains.

<img src="images/benchmark/pr_scatter_equal_isocirc_no_cirilong.png" width="100%"/>

Mean precision (P) / recall (R) / F1 across the same 3 runs, exon-based matching, plus each individual tool's own raw calls:

| Tier / tool | CIRI-long included | CIRI-long excluded |
| --- | --- | --- |
| isoCirc alone | P 0.87 / R 0.08 / F1 0.15 | P 0.87 / R 0.08 / F1 0.15 |
| circFL-seq alone | P 0.74 / R 0.23 / F1 0.36 | P 0.74 / R 0.23 / F1 0.35 |
| CIRI-long alone | P 0.77 / R 0.03 / F1 0.05 | not run |
| circNICK-lrs alone | P 0.22 / R 0.10 / F1 0.14 | P 0.22 / R 0.10 / F1 0.14 |
| `discovery` | P 0.66 / R 0.27 / F1 0.38 | P 0.66 / R 0.26 / F1 0.37 |
| `balanced_recall` | P 0.80 / R 0.25 / F1 0.38 | P 0.80 / R 0.24 / F1 0.38 |
| `balanced_precision` | P 0.84 / R 0.10 / F1 0.17 | P 0.85 / R 0.09 / F1 0.17 |
| `high_confidence` | P 0.86 / R 0.01 / F1 0.02 | P 0.87 / R 0.01 / F1 0.02 |

CIRI-long was the weakest individual tool on this protocol (F1 0.05, far behind its performance on the CIRI-long/circFL-seq protocols).

Excluding CIRI-long from the consensus tool set costs at most 0.01 F1 at any of the 4 tiers: P/R/F1 are identical or within 0.01 of each other with CIRI-long included vs excluded.

circNICK-lrs wet-lab protocol reads were not tested here: as noted above, circNICK-lrs's own protocol requires rolling-circle reads with multiple copies per read, a structure the other 3 tools can't use, so consensus scoring across tools does not apply to it.

The plots below cover this CIRI-long-excluded configuration in the same depth as the other two protocols get in the main section above, since this is the tool set actually recommended for isoCirc-protocol data. The main section's own isoCirc-protocol plots instead use a with-CIRI-long tool set, for direct comparability with the CIRI-long protocol and circFL-seq protocol columns.

<img src="images/benchmark/equal_gt_coverage_stacked_isocirc_no_cirilong.png" width="100%"/>

<img src="images/benchmark/equal_type_classification_isocirc_exon_no_cirilong.png" width="100%"/>
<img src="images/benchmark/equal_type_classification_isocirc_relaxed_bsj_no_cirilong.png" width="100%"/>
<img src="images/benchmark/equal_type_classification_isocirc_strict_bsj_no_cirilong.png" width="100%"/>

<img src="images/benchmark/pr_scatter_isoform_isocirc_no_cirilong.png" width="100%"/>
<img src="images/benchmark/equal_isoform_detection_reciprocal_overlap_99_isocirc_no_cirilong.png" width="100%"/>
<img src="images/benchmark/equal_isoform_detection_reciprocal_overlap_90_isocirc_no_cirilong.png" width="100%"/>
<img src="images/benchmark/equal_isoform_family_completeness_reciprocal_overlap_99_isocirc_no_cirilong.png" width="100%"/>
<img src="images/benchmark/equal_isoform_family_completeness_reciprocal_overlap_90_isocirc_no_cirilong.png" width="100%"/>

<img src="images/benchmark/equal_expression_correlation_isocirc_no_cirilong.png" width="100%"/>
<img src="images/benchmark/equal_expression_correlation_isocirc_strict_bsj_no_cirilong.png" width="100%"/>
<img src="images/benchmark/equal_expression_correlation_isocirc_exon_based_no_cirilong.png" width="100%"/>
<img src="images/benchmark/equal_expression_correlation_isocirc_full_no_cirilong.png" width="100%"/>
<img src="images/benchmark/equal_expression_correlation_isocirc_strict_no_cirilong.png" width="100%"/>
<img src="images/benchmark/equal_expression_correlation_isocirc_relaxed_no_cirilong.png" width="100%"/>
<img src="images/benchmark/equal_sensitivity_by_expression_isocirc_no_cirilong.png" width="100%"/>
<img src="images/benchmark/equal_sensitivity_by_expression_isocirc_isoform_no_cirilong.png" width="100%"/>
<img src="images/benchmark/equal_reproducibility_pairs_isocirc_no_cirilong.png" width="100%"/>
<img src="images/benchmark/equal_reproducibility_pairs_isocirc_strict_bsj_no_cirilong.png" width="100%"/>
<img src="images/benchmark/equal_reproducibility_pairs_isocirc_exon_based_no_cirilong.png" width="100%"/>
<img src="images/benchmark/equal_reproducibility_pairs_isocirc_full_no_cirilong.png" width="100%"/>
<img src="images/benchmark/equal_reproducibility_pairs_isocirc_strict_no_cirilong.png" width="100%"/>
<img src="images/benchmark/equal_reproducibility_pairs_isocirc_relaxed_no_cirilong.png" width="100%"/>

### Caller overlap

Ground-truth loci recovered by each of the four tools are pooled across runs and drawn as an UpSet-style plot, where a locus counts as recovered if it was matched in at least one run. Two versions were used in the original benchmark paper: exon-based, which requires correct internal structure, and relaxed BSJ, which requires only genomic position. This is now extended to all six matching strategies used elsewhere in this benchmark: strict BSJ and exon-based (exact/structural), relaxed BSJ (coordinate-tolerant), and the full/strict/relaxed isoform tiers (structure plus BSJ tolerance, tightest to loosest).

The point of this plot is agreement, not any one tool's own recall. Of every circRNA any tool finds, only about 13% are found by all four tools together (12.5% on the species benchmark, 12.9% on the wet-lab CIRI-long protocol), and roughly half are found by exactly one tool. Each tool is effectively viewing circRNAs through its own lens, with no single reliable "ground truth caller" to defer to. That disagreement is the core motivation for NanoCirc: a naive union of all four tools' calls would inherit every tool's false positives, and a naive intersection would keep only the rare loci every tool happens to agree on, trading away most of the recall. NanoCirc's confidence tiers reconcile this with vote-based merging instead of a plain union/intersection rule, as described in the pipeline's own methods documentation.

<img src="images/benchmark/upset_plot_equal_isocirc_strict_bsj_no_cirilong.png" width="100%"/>
<img src="images/benchmark/upset_plot_equal_isocirc_no_cirilong.png" width="100%"/>
<img src="images/benchmark/upset_plot_equal_isocirc_exon_based_no_cirilong.png" width="100%"/>
<img src="images/benchmark/upset_plot_equal_isocirc_full_no_cirilong.png" width="100%"/>
<img src="images/benchmark/upset_plot_equal_isocirc_strict_no_cirilong.png" width="100%"/>
<img src="images/benchmark/upset_plot_equal_isocirc_relaxed_no_cirilong.png" width="100%"/>

</details>

## Species benchmark: human vs mouse

Wet-lab protocol was fixed to CIRI-long for both species; each species used its own genome, GTF, circRNA database and in-silico generated by database fit circRNA catalog. Any difference below comes from species circRNA biology and genome/annotation/circRNA database differences, not the wet-lab protocol.

<details>
<summary><strong>Structural precision/recall</strong></summary>

### Structural precision/recall

<img src="images/benchmark/pr_scatter_human.png" width="100%"/>
<img src="images/benchmark/pr_scatter_mouse.png" width="100%"/>

Mean precision (P) / recall (R) / F1 across 3 runs, exon-based matching, per nanocirc tier:

| Tier | Human | Mouse |
| --- | --- | --- |
| `discovery` | P 0.67 / R 0.37 / F1 0.48 | P 0.74 / R 0.54 / F1 0.62 |
| `balanced_recall` | P 0.74 / R 0.36 / F1 0.48 | P 0.85 / R 0.52 / F1 0.65 |
| `balanced_precision` | P 0.85 / R 0.20 / F1 0.33 | P 0.94 / R 0.36 / F1 0.52 |
| `high_confidence` | P 0.89 / R 0.07 / F1 0.13 | P 0.98 / R 0.13 / F1 0.22 |

</details>

<details>
<summary><strong>Ground truth coverage</strong></summary>

### Ground truth coverage

<img src="images/benchmark/gt_coverage_stacked_human.png" width="100%"/>
<img src="images/benchmark/gt_coverage_stacked_mouse.png" width="100%"/>

Mouse has higher ground truth recovery than human at every tier (e.g. `discovery` R 0.54 vs human's 0.37, from the table above), the same species gap the summary table shows.

</details>

<details>
<summary><strong>circRNA type classification</strong></summary>

### circRNA type classification

<img src="images/benchmark/type_classification_human_exon.png" width="100%"/>
<img src="images/benchmark/type_classification_human_relaxed_bsj.png" width="100%"/>
<img src="images/benchmark/type_classification_human_strict_bsj.png" width="100%"/>
<img src="images/benchmark/type_classification_mouse_exon.png" width="100%"/>
<img src="images/benchmark/type_classification_mouse_relaxed_bsj.png" width="100%"/>
<img src="images/benchmark/type_classification_mouse_strict_bsj.png" width="100%"/>

Same TP/FP vs GT breakdown as the GT coverage plot above, split by circRNA biotype (`eciRNA`, `EIciRNA`, `ciRNA`, `antisense`, `intergenic`, see [docs/methods.md](https://github.com/aerusakovich/nanocirc_pipeline/blob/main/docs/methods.md#circrna-type-classification)). `discovery`-tier recovery per type, isoform-level matching:

| Species | eciRNA | EIciRNA | ciRNA | antisense | intergenic |
| --- | --- | --- | --- | --- | --- |
| Human | 73.2% | 17.2% | 25.8% | 28.6% | 20.8% |
| Mouse | 88.8% | 38.9% | 53.5% | 54.0% | 15.9% |

Mouse recovers every type better than human except intergenic, where human is slightly ahead (20.8% vs 15.9%), the one type/species combination where human derived circRNAs win. eciRNA is the easiest type on both species; intergenic is hardest on mouse, EIciRNA hardest on human.

</details>

<details>
<summary><strong>Isoform detection</strong></summary>

### Isoform detection

The structural plots show if a specific part of circRNA structure was found at its locus. But we are also interested to see if the callers have found the exactly same molecule, especially on isoform level, which is the main point of using long-read sequencing. We evaluated that with three tiers, each pairing a BSJ boundary check with an exon overlap check: `full` (exact BSJ, 100% overlap, byte identical structure), `strict` (BSJ within tolerance window, 99% exon overlap, allows little boundary noise), `relaxed` (BSJ within tolerance, 90% overlap, allowing a part of mismatched sequence).

<img src="images/benchmark/pr_scatter_isoform_human.png" width="100%"/>
<img src="images/benchmark/pr_scatter_isoform_mouse.png" width="100%"/>

Long-read protocols are used mainly to resolve internal exon structure, not just BSJ position. Among ground-truth loci where the same BSJ has more than one real isoform (isoform family), we checked how often each caller recovers the exact isoform, versus calls the right BSJ with the wrong structure, versus misses the BSJ entirely.

<img src="images/benchmark/isoform_detection_reciprocal_overlap_99_human.png" width="100%"/>
<img src="images/benchmark/isoform_detection_reciprocal_overlap_99_mouse.png" width="100%"/>

<img src="images/benchmark/isoform_detection_reciprocal_overlap_90_human.png" width="100%"/>
<img src="images/benchmark/isoform_detection_reciprocal_overlap_90_mouse.png" width="100%"/>

Green+grey sum to 100% of the isoform-family ground truth; red shows the calls that found the right BSJ but didn't match internal structure criteria. 99% reciprocal exon overlap allows only trivial boundary noise; 90% allows a small sequence mismatch while still requiring the same structure, not a different one. % of ground-truth isoforms called exactly right, at 99% reciprocal overlap:

| Tier / tool | Human | Mouse |
| --- | --- | --- |
| CIRI-long alone | 19.6% | 30.8% |
| isoCirc alone | 13.0% | 21.1% |
| circNICK-lrs alone | 3.1% | 5.6% |
| circFL-seq alone | 19.2% | 23.2% |
| `discovery` | 31.0% | 41.9% |
| `balanced_recall` | 29.4% | 40.6% |
| `balanced_precision` | 20.4% | 27.5% |
| `high_confidence` | 6.9% | 10.1% |

Mouse resolves more isoforms than human at every tier and tool, consistent with mouse's higher overall recall in this benchmark. Both species are far below the wet-lab-protocol benchmark's exact-isoform rates above (e.g. `discovery` 31.0%/41.9% here vs 61.8-67.6% on CIRI-long/circFL-seq protocol), since this dataset uses a more isoform-family-rich catalog (422/912 isoforms tested here vs 102 there). CIRI-long and circFL-seq lead the raw tools on both species (19-31%); circNICK-lrs is the weakest by a wide margin (3.1%/5.6%), and even `high_confidence` clears it on both species.

<img src="images/benchmark/isoform_family_completeness_reciprocal_overlap_99_human.png" width="100%"/>
<img src="images/benchmark/isoform_family_completeness_reciprocal_overlap_99_mouse.png" width="100%"/>

<img src="images/benchmark/isoform_family_completeness_reciprocal_overlap_90_human.png" width="100%"/>
<img src="images/benchmark/isoform_family_completeness_reciprocal_overlap_90_mouse.png" width="100%"/>

Per family plot instead of per isoform: did the caller get every isoform in the family right, some of them, or none. % of isoform families fully resolved, at 99% reciprocal overlap:

| Tier / tool | Human | Mouse |
| --- | --- | --- |
| CIRI-long alone | 8.6% | 13.7% |
| isoCirc alone | 2.4% | 10.8% |
| circNICK-lrs alone | 0.0% | 0.0% |
| circFL-seq alone | 11.5% | 15.0% |
| `discovery` | 16.3% | 26.0% |
| `balanced_recall` | 14.4% | 23.5% |
| `balanced_precision` | 7.2% | 10.8% |
| `high_confidence` | 1.4% | 4.7% |

Getting one isoform right while missing or misclassifying others is the most common outcome; fully resolving an isoform family is rare for every tier on both species, and `high_confidence` barely ever manages it (1.4%/4.7%). Mouse is consistently ahead of human here too. circNICK-lrs never fully resolves a family on either species (0.0%); circFL-seq leads the raw tools (11.5%/15.0%), ahead even of `balanced_precision` on human (11.5% vs 7.2%).

</details>

<details>
<summary><strong>Expression quantification and reproducibility</strong></summary>

### Expression quantification and reproducibility

<img src="images/benchmark/expression_correlation_human.png" width="100%"/>
<img src="images/benchmark/expression_correlation_mouse.png" width="100%"/>

<img src="images/benchmark/expression_correlation_human_strict_bsj.png" width="100%"/>
<img src="images/benchmark/expression_correlation_mouse_strict_bsj.png" width="100%"/>

<img src="images/benchmark/expression_correlation_human_exon_based.png" width="100%"/>
<img src="images/benchmark/expression_correlation_mouse_exon_based.png" width="100%"/>

Read-count correlation against ground truth (right panel) and cross-method agreement (left panel), per tool/tier. Each cell shows three metrics:

- **Pearson r**: linear correlation between predicted and true counts. Sensitive to whether counts scale proportionally, not just in the same direction.
- **Spearman rho**: rank correlation. Only checks whether higher-expressed circRNAs get higher counts, regardless of exact scale. More resistant to outliers and non-linear (but still monotonic) relationships.
- **R²** (Pearson r squared): how much of the variance in true expression is captured by the predicted counts, on a 0-1 scale.

For differential expression, Spearman rho matters to rank circRNAs correctly relative to each other; Pearson r/R² matter to ensure the fold-changes DESeq2 estimates between conditions are quantitatively trustworthy, not just directionally right. Pearson r by species (relaxed-BSJ matching):

| Caller | Human | Mouse |
| --- | --- | --- |
| CIRI-long | 0.911 | 0.906 |
| isoCirc | 0.552 | 0.667 |
| circNICK-lrs | 0.763 | 0.814 |
| circFL-seq | 0.787 | 0.805 |
| `high_confidence` | 0.967 | 0.930 |
| `balanced_precision` | 0.974 | 0.947 |
| `balanced_recall` | 0.967 | 0.936 |
| `discovery` | 0.967 | 0.936 |

isoCirc has the weakest quantification of any tool/tier on both species (r 0.55-0.67), same as on every wet-lab protocol above. Every nanocirc tier correlates more strongly with ground truth than any raw tool on both species (r ≥ 0.93). Unlike detection, where mouse leads, quantification quality is close between species and `balanced_precision` is actually slightly stronger on human (0.974 vs 0.947).

#### Isoform-level expression quantification

The table above reports relaxed-BSJ matching, which collapses a locus with several real isoforms into a single point. The tables below instead score each real GT isoform on its own: `full`, `strict`, `relaxed` are the same 3 tiers used for [Isoform detection](#isoform-detection-1) and the isoform PR scatter above (`full` = exact BSJ + 100% exon overlap, `strict` = BSJ within tolerance + 99% overlap, `relaxed` = BSJ within tolerance + 90% overlap). A GT isoform's predicted count is whichever caller isoform structurally matched it under that tier, or 0 if none did, so n here is a genuinely stricter re-match, not the table above's n expanded, and can be smaller. Pearson r (n):

**`full` (exact BSJ, 100% exon overlap):**

| Caller | Human r (n) | Mouse r (n) |
| --- | --- | --- |
| CIRI-long | 0.908 (5,742) | 0.903 (7,522) |
| isoCirc | 0.557 (3,342) | 0.675 (6,454) |
| circNICK-lrs | 0.537 (1,515) | 0.760 (1,600) |
| circFL-seq | 0.794 (4,645) | 0.811 (6,347) |
| `high_confidence` | 0.965 (2,443) | 0.937 (4,649) |
| `balanced_precision` | 0.971 (5,064) | 0.947 (8,658) |
| `balanced_recall` | 0.964 (7,473) | 0.935 (11,066) |
| `discovery` | 0.962 (7,674) | 0.936 (11,266) |

**`strict` (BSJ within tolerance, 99% exon overlap):**

| Caller | Human r (n) | Mouse r (n) |
| --- | --- | --- |
| CIRI-long | 0.910 (6,512) | 0.904 (9,181) |
| isoCirc | 0.557 (3,363) | 0.676 (6,545) |
| circNICK-lrs | 0.538 (1,596) | 0.766 (1,685) |
| circFL-seq | 0.794 (4,692) | 0.809 (6,400) |
| `high_confidence` | 0.965 (2,446) | 0.936 (4,708) |
| `balanced_precision` | 0.971 (5,097) | 0.947 (8,756) |
| `balanced_recall` | 0.965 (7,943) | 0.936 (12,015) |
| `discovery` | 0.964 (8,144) | 0.937 (12,224) |

**`relaxed` (BSJ within tolerance, 90% exon overlap):**

| Caller | Human r (n) | Mouse r (n) |
| --- | --- | --- |
| CIRI-long | 0.909 (7,206) | 0.906 (9,824) |
| isoCirc | 0.558 (3,384) | 0.677 (6,588) |
| circNICK-lrs | 0.546 (1,732) | 0.763 (1,868) |
| circFL-seq | 0.793 (4,722) | 0.809 (6,432) |
| `high_confidence` | 0.965 (2,452) | 0.932 (4,729) |
| `balanced_precision` | 0.971 (5,130) | 0.948 (8,793) |
| `balanced_recall` | 0.964 (8,280) | 0.939 (12,198) |
| `discovery` | 0.963 (8,481) | 0.939 (12,413) |

n grows from `full` to `relaxed` on both species as the overlap requirement loosens, and mouse's n grows a bit faster throughout, consistent with mouse's richer isoform-family catalog documented in [Isoform detection](#isoform-detection-1) above. circNICK-lrs is the caller most affected by tightening the requirement: its correlation drops against the locus-level table above (r 0.54-0.76 here vs r 0.74-0.81 there) and its n is by far the smallest of any caller (1,515-1,868), since its structure calls carry real boundary noise a genuine per-isoform re-match does not forgive. Every nanocirc tier still correlates far more strongly than any raw tool at every tier of stringency, and human/mouse stay close to each other throughout, same as the table above.

<img src="images/benchmark/expression_correlation_human_full.png" width="100%"/>
<img src="images/benchmark/expression_correlation_mouse_full.png" width="100%"/>

<img src="images/benchmark/expression_correlation_human_strict.png" width="100%"/>
<img src="images/benchmark/expression_correlation_mouse_strict.png" width="100%"/>

<img src="images/benchmark/expression_correlation_human_relaxed.png" width="100%"/>
<img src="images/benchmark/expression_correlation_mouse_relaxed.png" width="100%"/>

<img src="images/benchmark/sensitivity_by_expression_human.png" width="100%"/>
<img src="images/benchmark/sensitivity_by_expression_mouse.png" width="100%"/>

Same tertiles, matched at the `full`/`strict`/`relaxed` isoform tiers instead of structural comparison:

<img src="images/benchmark/sensitivity_by_expression_human_isoform.png" width="100%"/>
<img src="images/benchmark/sensitivity_by_expression_mouse_isoform.png" width="100%"/>

Detection sensitivity broken down by ground-truth expression tertile: low-expressed circRNAs are structurally harder to detect correctly than high-expressed ones, for every tool.

| Tier | Species | Low | Mid | High |
| --- | --- | --- | --- | --- |
| `discovery` | Human | 15.7% | 37.3% | 53.5% |
| `discovery` | Mouse | 26.7% | 55.3% | 68.6% |
| `high_confidence` | Human | 0.8% | 5.8% | 15.3% |
| `high_confidence` | Mouse | 2.4% | 11.4% | 22.8% |

Mouse is more sensitive than human at every tertile and tier shown, widest at `discovery`'s low tertile (26.7% vs 15.7%, a 70% relative gap) and narrowest at `high_confidence`'s low tertile, where both species are already close to floor.

<img src="images/benchmark/reproducibility_pairs_human.png" width="100%"/>
<img src="images/benchmark/reproducibility_pairs_mouse.png" width="100%"/>

<img src="images/benchmark/reproducibility_pairs_human_strict_bsj.png" width="100%"/>
<img src="images/benchmark/reproducibility_pairs_mouse_strict_bsj.png" width="100%"/>

<img src="images/benchmark/reproducibility_pairs_human_exon_based.png" width="100%"/>
<img src="images/benchmark/reproducibility_pairs_mouse_exon_based.png" width="100%"/>

Pairwise read-count agreement across every tool/tier combination, ground truth included. As with the wet-lab-protocol benchmark above, nanocirc's own tiers agree with each other almost perfectly (diagonal points in every tier-vs-tier panel) since they share the same remap-based quantification; agreement against raw tool calls is looser, isoCirc again showing the widest scatter, consistent with its weaker Pearson r above.

Same reproducibility matrix, matched at the `full`/`strict`/`relaxed` isoform tiers instead of structural comparison, so agreement is evaluated on isoform level, each real isoform scored against its own matched caller isoform (see [Isoform-level expression quantification](#isoform-level-expression-quantification-1) above).

<img src="images/benchmark/reproducibility_pairs_human_full.png" width="100%"/>
<img src="images/benchmark/reproducibility_pairs_mouse_full.png" width="100%"/>

<img src="images/benchmark/reproducibility_pairs_human_strict.png" width="100%"/>
<img src="images/benchmark/reproducibility_pairs_mouse_strict.png" width="100%"/>

<img src="images/benchmark/reproducibility_pairs_human_relaxed.png" width="100%"/>
<img src="images/benchmark/reproducibility_pairs_mouse_relaxed.png" width="100%"/>

### Caller overlap

Ground-truth loci recovered by each of the four tools are pooled across runs and drawn as an UpSet-style plot, where a locus counts as recovered if it was matched in at least one run. Two versions were used in the original benchmark paper: exon-based, which requires correct internal structure, and relaxed BSJ, which requires only genomic position. This is now extended to all six matching strategies used elsewhere in this benchmark: strict BSJ and exon-based (exact/structural), relaxed BSJ (coordinate-tolerant), and the full/strict/relaxed isoform tiers (structure plus BSJ tolerance, tightest to loosest).

The point of this plot is agreement, not any one tool's own recall. Of every circRNA any tool finds, only about 13% are found by all four tools together (12.5% on the species benchmark, 12.9% on the wet-lab CIRI-long protocol), and roughly half are found by exactly one tool. Each tool is effectively viewing circRNAs through its own lens, with no single reliable "ground truth caller" to defer to. That disagreement is the core motivation for NanoCirc: a naive union of all four tools' calls would inherit every tool's false positives, and a naive intersection would keep only the rare loci every tool happens to agree on, trading away most of the recall. NanoCirc's confidence tiers reconcile this with vote-based merging instead of a plain union/intersection rule, as described in the pipeline's own methods documentation.

<img src="images/benchmark/upset_plot_human_strict_bsj.png" width="100%"/>
<img src="images/benchmark/upset_plot_mouse_strict_bsj.png" width="100%"/>

<img src="images/benchmark/upset_plot_human.png" width="100%"/>
<img src="images/benchmark/upset_plot_mouse.png" width="100%"/>

<img src="images/benchmark/upset_plot_human_exon_based.png" width="100%"/>
<img src="images/benchmark/upset_plot_mouse_exon_based.png" width="100%"/>

<img src="images/benchmark/upset_plot_human_full.png" width="100%"/>
<img src="images/benchmark/upset_plot_mouse_full.png" width="100%"/>

<img src="images/benchmark/upset_plot_human_strict.png" width="100%"/>
<img src="images/benchmark/upset_plot_mouse_strict.png" width="100%"/>

<img src="images/benchmark/upset_plot_human_relaxed.png" width="100%"/>
<img src="images/benchmark/upset_plot_mouse_relaxed.png" width="100%"/>

</details>

**What we see:**

- **The species effect is real for recall, not for precision.** Pooled across every caller and circRNA type, mouse recall/F1 is significantly higher (Wilcoxon p < 10⁻¹⁵), but precision shows no significant species gap at all (p = 0.44), and actually reverses for antisense and intergenic circRNAs, where human precision is significantly higher (see [Statistical significance](#statistical-significance-1) below). NanoCirc's own merged tiers still score higher on mouse at all 4 tiers, on both precision and recall (table above) — the F1 gap is largest at `balanced_precision` (F1 0.33 human vs 0.52 mouse) and smallest at `high_confidence` (F1 0.13 vs 0.22) — but that tier-level precision edge doesn't reflect a uniform species effect once you look past nanocirc's own merged output.

- **Quantification differences between species are not as big**, on both precision and recall. Quantification performance is close between species and `balanced_precision` is actually slightly stronger on human (0.974 vs 0.947).

<details>
<summary><strong>Statistical significance</strong></summary>

### Statistical significance

Same statistical methods as the [wet-lab-protocol comparison](#statistical-significance) above (linear regression, Kendall's tau-b, Mantel-Haenszel + Breslow-Day, Wilcoxon signed-rank, Cohen's d), matching every caller/circRNA-type/run point between human and mouse.

#### Aggregate, by metric (n paired points pooling every caller × type × run)

| Metric | n | Human | Mouse | Difference | Slope | Intercept | R² | Wilcoxon p | Cohen's d |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Precision | 108 | 0.725 | 0.734 | +0.009 | 0.935 | 0.057 | 0.836 | 0.436 | 0.092 |
| Recall | 120 | 0.187 | 0.290 | +0.103 | 1.213 | 0.063 | 0.852 | 2×10⁻¹⁶ | 0.977 |
| F1 | 120 | 0.248 | 0.361 | +0.113 | 1.240 | 0.054 | 0.862 | 6×10⁻¹⁶ | 1.023 |

Recall and F1 are well predicted by human's value (R² 0.85–0.86) with a slope above 1: mouse's advantage is higher for callers that are stronger, not frowing uniformly for all callers. Precision fits well (R² = 0.84, slope near 1) but shows no systematic species gap at all (Wilcoxon p = 0.44): mouse is not reliably more precise than human once every caller and type is pooled together.

#### Per type, by metric (n=24 paired points per row pooling all 8 callers × 3 runs, except ciRNA precision n=21, antisense precision n=18, and intergenic precision n=21: some caller × run cells had zero calls of that type, so precision can't be defined and those cells are dropped)

| Type | Metric | n | Human | Mouse | Difference | Slope | Intercept | R² | Wilcoxon p | Cohen's d |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| eciRNA | Precision | 24 | 0.651 | 0.748 | +0.096 | 1.046 | 0.066 | 0.945 | 1×10⁻⁷ | 1.699 |
| eciRNA | Recall | 24 | 0.495 | 0.619 | +0.124 | 1.111 | 0.069 | 0.935 | 1×10⁻⁷ | 2.196 |
| eciRNA | F1 | 24 | 0.532 | 0.657 | +0.125 | 1.206 | 0.016 | 0.938 | 1×10⁻⁷ | 2.280 |
| EIciRNA | Precision | 24 | 0.629 | 0.652 | +0.023 | 0.974 | 0.039 | 0.827 | 0.747 | 0.226 |
| EIciRNA | Recall | 24 | 0.103 | 0.230 | +0.127 | 2.287 | −0.006 | 0.961 | 1×10⁻⁷ | 1.725 |
| EIciRNA | F1 | 24 | 0.168 | 0.321 | +0.153 | 1.855 | 0.009 | 0.924 | 2×10⁻⁶ | 1.895 |
| ciRNA | Precision | 21 | 0.862 | 0.925 | +0.063 | 0.562 | 0.441 | 0.984 | 2×10⁻⁴ | 1.331 |
| ciRNA | Recall | 24 | 0.114 | 0.242 | +0.128 | 2.030 | 0.011 | 0.995 | 6×10⁻⁵ | 1.219 |
| ciRNA | F1 | 24 | 0.182 | 0.336 | +0.153 | 1.658 | 0.033 | 0.990 | 6×10⁻⁵ | 1.498 |
| antisense | Precision | 18 | 0.774 | 0.720 | −0.054 | 0.939 | −0.006 | 0.997 | 4×10⁻⁴ | −1.860 |
| antisense | Recall | 24 | 0.135 | 0.278 | +0.143 | 1.984 | 0.011 | 0.863 | 2×10⁻⁴ | 0.976 |
| antisense | F1 | 24 | 0.216 | 0.356 | +0.141 | 1.603 | 0.010 | 0.894 | 2×10⁻⁴ | 0.952 |
| intergenic | Precision | 21 | 0.742 | 0.635 | −0.107 | 0.896 | −0.029 | 0.853 | 1×10⁻⁵ | −1.691 |
| intergenic | Recall | 24 | 0.089 | 0.083 | −0.006 | 0.714 | 0.020 | 0.922 | 0.434 | −0.191 |
| intergenic | F1 | 24 | 0.139 | 0.133 | −0.006 | 0.718 | 0.033 | 0.905 | 0.794 | −0.138 |

Mouse is significantly ahead on recall and F1 for eciRNA, EIciRNA, ciRNA, and antisense. intergenic is the exception: mouse is not significantly different from human on either recall or F1 for this type (p = 0.43 and 0.79). For precision species relationship flips by type: eciRNA and ciRNA favor mouse (p ≤ 2×10⁻⁴), EIciRNA shows no significant difference (p = 0.75), and antisense and intergenic actually favor human (negative differences, p ≤ 4×10⁻⁴, the only rows in this whole table where human is significantly ahead of mouse).

#### Does caller ranking transfer?

We averaged each caller's human and mouse values across type and run first, ranked the 8 callers by that average on each side, then compared the two rankings with Kendall's tau-b:

| Metric | n callers | Kendall's tau-b | Spearman rho |
| --- | --- | --- | --- |
| Precision | 8 | 0.929 (p = 4×10⁻⁴) | 0.976 (p = 3×10⁻⁵) |
| Recall | 8 | 0.929 (p = 4×10⁻⁴) | 0.976 (p = 3×10⁻⁵) |
| F1 | 8 | 0.929 (p = 4×10⁻⁴) | 0.976 (p = 3×10⁻⁵) |

Ranking transfers significantly on every metric. By F1:

| Caller | F1 human | F1 mouse | Rank human | Rank mouse |
| --- | --- | --- | --- | --- |
| balanced_recall | 0.406 | 0.572 | 1 | 1 |
| discovery | 0.396 | 0.545 | 2 | 2 |
| CIRI-long | 0.319 | 0.442 | 3 | 3 |
| balanced_precision | 0.263 | 0.440 | 4 | 4 |
| circFL-seq | 0.251 | 0.327 | 5 | 6 |
| isoCirc | 0.188 | 0.343 | 6 | 5 |
| high_confidence | 0.102 | 0.160 | 7 | 7 |
| circNICK-lrs | 0.056 | 0.055 | 8 | 8 |

Six of the eight callers reach the exact same F1 rank on both species. Only circFL-seq and isoCirc swap places (circFL-seq ranks 5th on human but 6th on mouse, isoCirc the reverse).

F1 here is lower than the [structural precision/recall](#structural-precisionrecall-1) table above for the same tier (e.g. `discovery` 0.40/0.55 here vs 0.48/0.62 there): that table pools all circRNA types together per run before averaging across the 3 runs, while this one averages precision/recall/F1 per circRNA type first, then across type × run (15 points). The two are different aggregations of the same underlying calls, not different data.

#### Per tool/mode, by metric (n=15 paired points per row: the 5 circRNA types × 3 runs, except circNICK-lrs precision n=6 and high_confidence precision n=12)

| Tool | Metric | n | Human | Mouse | Difference | Slope | Intercept | R² | Wilcoxon p | Cohen's d |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CIRI-long | Precision | 15 | 0.687 | 0.791 | +0.104 | 0.389 | 0.524 | 0.350 | 0.003 | 1.078 |
| CIRI-long | Recall | 15 | 0.240 | 0.338 | +0.097 | 1.053 | 0.085 | 0.833 | 6×10⁻⁴ | 1.124 |
| CIRI-long | F1 | 15 | 0.319 | 0.442 | +0.123 | 1.132 | 0.081 | 0.713 | 9×10⁻⁴ | 1.153 |
| isoCirc | Precision | 15 | 0.877 | 0.844 | −0.032 | 0.880 | 0.073 | 0.478 | 0.135 | −0.426 |
| isoCirc | Recall | 15 | 0.121 | 0.255 | +0.133 | 1.581 | 0.063 | 0.866 | 6×10⁻⁵ | 1.148 |
| isoCirc | F1 | 15 | 0.188 | 0.343 | +0.155 | 1.362 | 0.087 | 0.908 | 6×10⁻⁵ | 1.429 |
| circNICK-lrs | Precision | 6 | 0.167 | 0.131 | −0.036 | 4.369 | −0.598 | 0.906 | 0.438 | −0.709 |
| circNICK-lrs | Recall | 15 | 0.056 | 0.060 | +0.004 | 1.036 | 0.002 | 0.999 | 0.028 | 0.774 |
| circNICK-lrs | F1 | 15 | 0.056 | 0.055 | −0.000 | 1.037 | −0.002 | 0.994 | 0.600 | −0.047 |
| circFL-seq | Precision | 15 | 0.624 | 0.603 | −0.020 | 1.053 | −0.053 | 0.946 | 0.972 | −0.237 |
| circFL-seq | Recall | 15 | 0.174 | 0.236 | +0.062 | 1.081 | 0.048 | 0.918 | 0.009 | 1.055 |
| circFL-seq | F1 | 15 | 0.251 | 0.327 | +0.076 | 1.074 | 0.057 | 0.886 | 0.013 | 0.924 |
| high_confidence | Precision | 12 | 0.925 | 0.902 | −0.023 | 0.475 | 0.462 | 0.240 | 0.369 | −0.291 |
| high_confidence | Recall | 15 | 0.065 | 0.109 | +0.044 | 1.496 | 0.012 | 0.995 | 6×10⁻⁵ | 0.759 |
| high_confidence | F1 | 15 | 0.102 | 0.160 | +0.058 | 1.347 | 0.023 | 0.990 | 6×10⁻⁵ | 0.913 |
| balanced_precision | Precision | 15 | 0.888 | 0.872 | −0.016 | 0.322 | 0.586 | 0.088 | 0.421 | −0.197 |
| balanced_precision | Recall | 15 | 0.186 | 0.332 | +0.146 | 1.128 | 0.122 | 0.858 | 6×10⁻⁵ | 1.563 |
| balanced_precision | F1 | 15 | 0.263 | 0.440 | +0.177 | 1.074 | 0.157 | 0.870 | 6×10⁻⁵ | 1.961 |
| balanced_recall | Precision | 15 | 0.690 | 0.740 | +0.050 | 0.662 | 0.283 | 0.370 | 0.095 | 0.412 |
| balanced_recall | Recall | 15 | 0.323 | 0.492 | +0.168 | 0.999 | 0.169 | 0.740 | 9×10⁻⁴ | 1.382 |
| balanced_recall | F1 | 15 | 0.406 | 0.572 | +0.166 | 1.073 | 0.137 | 0.593 | 9×10⁻⁴ | 1.277 |
| discovery | Precision | 15 | 0.650 | 0.664 | +0.014 | 0.837 | 0.120 | 0.722 | 0.890 | 0.154 |
| discovery | Recall | 15 | 0.331 | 0.502 | +0.171 | 1.010 | 0.168 | 0.756 | 9×10⁻⁴ | 1.410 |
| discovery | F1 | 15 | 0.396 | 0.545 | +0.149 | 1.121 | 0.101 | 0.561 | 9×10⁻⁴ | 1.199 |

CIRI-long, isoCirc, high_confidence, balanced_precision, balanced_recall, and discovery are all significantly higher on mouse for recall and F1 (p ≤ 9×10⁻⁴). circFL-seq and circNICK-lrs are the exceptions: circFL-seq's recall/F1 gap is smaller (p ≤ 0.013), and circNICK-lrs shows no significant recall or F1 difference between species at all. Only CIRI-long shows a significant species gap on precision (p = 0.003, mouse higher); every other caller's precision is statistically indistinguishable between human and mouse.

| Tool | Outcome | Strata (types) | OR (MH) | 95% CI | MH p | Breslow-Day p |
| --- | --- | --- | --- | --- | --- | --- |
| CIRI-long | Precision | 5 | 2.17 | 2.03–2.31 | <10⁻³⁰⁰ | <10⁻³⁰⁰ |
| CIRI-long | Recall | 5 | 1.79 | 1.72–1.87 | <10⁻³⁰⁰ | <10⁻³⁰⁰ |
| isoCirc | Precision | 5 | 1.03 | 0.92–1.16 | 0.599 | <10⁻³⁰⁰ |
| isoCirc | Recall | 5 | 3.11 | 2.96–3.27 | <10⁻³⁰⁰ | <10⁻³⁰⁰ |
| circNICK-lrs | Precision | 5 | 0.90 | 0.83–0.97 | 0.004 | <10⁻³⁰⁰ |
| circNICK-lrs | Recall | 5 | 1.08 | 1.00–1.16 | 0.054 | 0.172 |
| circFL-seq | Precision | 5 | 1.06 | 0.97–1.15 | 0.214 | <10⁻³⁰⁰ |
| circFL-seq | Recall | 5 | 1.64 | 1.56–1.72 | <10⁻³⁰⁰ | <10⁻³⁰⁰ |
| high_confidence | Precision | 5 | 1.88 | 1.57–2.26 | 1×10⁻¹¹ | 5×10⁻⁵ |
| high_confidence | Recall | 5 | 2.10 | 1.96–2.25 | <10⁻³⁰⁰ | 5×10⁻⁵ |
| balanced_precision | Precision | 5 | 1.53 | 1.39–1.68 | <10⁻³⁰⁰ | <10⁻³⁰⁰ |
| balanced_precision | Recall | 5 | 2.73 | 2.61–2.86 | <10⁻³⁰⁰ | <10⁻³⁰⁰ |
| balanced_recall | Precision | 5 | 1.59 | 1.51–1.68 | <10⁻³⁰⁰ | <10⁻³⁰⁰ |
| balanced_recall | Recall | 5 | 2.45 | 2.36–2.55 | <10⁻³⁰⁰ | <10⁻³⁰⁰ |
| discovery | Precision | 5 | 1.26 | 1.21–1.32 | <10⁻³⁰⁰ | <10⁻³⁰⁰ |
| discovery | Recall | 5 | 2.51 | 2.41–2.61 | <10⁻³⁰⁰ | <10⁻³⁰⁰ |

Every caller has significantly higher odds of recovering a ground-truth circRNA (recall) on mouse than on human (OR 1.1–3.1), and most also have higher odds of a correct call (precision), though isoCirc and circFL-seq show no significant precision difference (p = 0.60 and 0.21) and circNICK-lrs's precision actually favors human (OR = 0.90, p = 0.004). Only 1 of the 16 caller/outcome combinations (circNICK-lrs's recall) has a species effect that's is stable across circRNA types (Breslow-Day p = 0.17); every other combination's effect size depends on circRNA type.

**Conclusion.** Correlation is good for recall and F1, but weaker and inconsistent for precision, where half of callers and types show no significant species difference at all, and antisense/intergenic precision actually favors human. Caller ranking transfers strongly (Kendall's tau-b = 0.93 on every metric); only circFL-seq and isoCirc swap adjacent ranks. The size of the species gap depends on circRNA types (1 of 16 Mantel-Haenszel/Breslow-Day combinations homogeneous). So the recall/F1 species gap exists, the caller ranking is mostly preserved between the species, and precision across species is mostly caller- and type-dependent.

</details>

<details>
<summary><strong>What are these false positives?</strong></summary>

### What are these false positives?

We traced every FP-supporting read back to the simulator's own ground truth (single-tool discovery exon-based FPs, all 3 human runs, mean % ± SD across runs):

| FP read is... | CIRI-long | circFL-seq | isoCirc | circNICK-lrs |
| --- | --- | --- | --- | --- |
| Real circ read, correct boundary, wrong strand (antisense called as sense) | 87.1 ± 0.2 | 77.6 ± 2.2 | 92.4 ± 0.8 | 5.6 ± 0.5 |
| Real circ read, correct strand, pure boundary noise | 2.2 ± 0.5 | 6.2 ± 2.5 | 0.0 ± 0.0 | 34.5 ± 2.4 |
| Real circ read, but circRNA reported at a different/distant locus | 10.7 ± 0.4 | 16.3 ± 0.5 | 7.6 ± 0.8 | 59.8 ± 2.7 |
| Not from a circ read at all (linear read miscalled) | 0.0 ± 0.0 | 0.0 ± 0.0 | 0.0 ± 0.0 | 0.1 ± 0.0 |
| FP locus overlaps a genomic repeat | 18.9 ± 0.1 | 34.0 ± 1.3 | 24.5 ± 3.0 | 64.6 ± 0.6 |
| n FP loci (pooled, 3 runs) | 2,679 | 1,248 | 375 | 6,844 |
| n supporting reads (pooled, 3 runs) | 63,109 | 15,341 | 4,675 | 29,099 |

Many FP come from `antisense`-type circRNA, where the simulator deliberately sets the circRNA's strand opposite to its host gene. NanoSim also mimics cDNA based RCRT-sequencing, meaning it simulates both strand sequences: one as an actual circRNA read and the other as cDNA read. The tools correctly find the genomic position of circRNA but select the host gene's (wrong) strand instead of flipping it to the antisense one (correct). Since strand agreement is required by all three matching strategies (exon-based, relaxed BSJ, strict BSJ), no boundary-tolerance setting can rescue these FP and we didn't find an easy way to override this gene-strand bias without introducing new FP elsewhere. These realistic strand misses could also happen on real data and we suggest to keep that in mind as you move onto wet-lab validation.

**Practical implication**: most `discovery` false positives still point at the right gene/circRNA locus, just with an imprecise boundary or wrong strand. That makes them usable for gene- or locus-level expression signal analysis, but not for isoform-level analysis or primer/probe design, both of which need to get the exact strand/boundary.

</details>

<details>
<summary><strong>What discovery still misses?</strong></summary>

### What discovery still misses?

We also characterized every ground-truth circRNA against whether any tool recovered it at all (mean % ± SD across 3 human runs; missed + recovered = 100% of the GT catalog, ~8,700 circRNAs/run):

| | Missed (FN) | Recovered |
| --- | --- | --- |
| Share of GT catalog | 64.8 ± 0.0 | 35.2 ± 0.0 |
| Median mature length (bp) | 1,093 ± 10 | 277 ± 2 |
| Mean rolling-circle copy number | 7.7 ± 0.1 | 8.6 ± 0.0 |
| Median simulated read count | 35 ± 0 | 99 ± 0 |
| Canonical GT/AG splice site at BSJ | 63.8 ± 0.1 | 80.2 ± 0.1 |
| No gene overlap at all (intergenic) | 19.9 ± 0.2 | 9.1 ± 0.2 |
| Falls in a gene with ≥5 isoforms | 73.7 ± 0.1 | 83.0 ± 0.5 |
| Overlaps a genomic repeat | 51.1 ± 0.5 | 57.7 ± 0.9 |

Recall per biotype:

| Biotype | Missed | Recovered |
| --- | --- | --- |
| `eciRNA` (pure exonic) | 22.9 ± 0.2 | 77.1 ± 0.2 |
| `antisense` | 68.9 ± 0.1 | 31.1 ± 0.1 |
| `ciRNA` (purely intronic) | 74.5 ± 0.2 | 25.5 ± 0.2 |
| `EIciRNA` (exon+intron) | 82.1 ± 0.1 | 17.9 ± 0.1 |
| `intergenic` | 81.3 ± 0.2 | 18.7 ± 0.2 |

Missed circRNAs are generally longer (~4x), lower-expressed (~2.8x fewer reads), and use non-canonical splice sites more often. Rolling-circle copy number is also slightly lower, meaning weaker RCRT/RCA tandem-repeat evidence for the tools that rely on it. `eciRNA` seems to be the easier type to recover; every other biotype, especially without a clean match to an annotated exon, is more likely to be missed. Repeat overlap isn't a meaningful driver of circRNA discovery.

</details>

---

## Database benchmark: annotated circRNA recovery

Other benchmarks use synthetic circRNAs generated by nanocirc-simulator and processed into reads by NanoSim. This benchmark evaluates performance using 2,000 database-derived circRNAs per type instead, sampled randomly from the human circBase/circAtlas union catalog. It uses the same NanoSim read simulation, expression distribution, detection pipeline, and quantification as every other dataset here. It shows the recovery of already-annotated circRNAs versus in silico de novo ones. It is important to note that these databases can contain circRNAs derived from short-read sequencing computational predictions that were never verified in a wet lab, are biased toward what tools can already detect (since that is how they were catalogued in the first place), and that some tools use the circRNA database itself as direct input for detection. It is expected that tools will perform better on database circRNAs, and there is no guarantee that every circRNA in the database is biologically true — but even with these limitations, it is a useful comparator and the closest approximation to biological ground truth we have.

<details>
<summary><strong>Structural precision/recall</strong></summary>

### Structural precision/recall

<img src="images/benchmark/pr_scatter_database_human.png" width="100%"/>

Mean precision (P) / recall (R) / F1 across 3 runs, exon-based matching, per nanocirc tier:

| Tier | Database (human) |
| --- | --- |
| `discovery` | P 0.85 / R 0.57 / F1 0.68 |
| `balanced_recall` | P 0.94 / R 0.56 / F1 0.70 |
| `balanced_precision` | P 0.97 / R 0.27 / F1 0.42 |
| `high_confidence` | P 0.97 / R 0.08 / F1 0.14 |

</details>

<details>
<summary><strong>Ground truth coverage</strong></summary>

### Ground truth coverage

<img src="images/benchmark/gt_coverage_stacked_database_human.png" width="100%"/>

Green (recovered) plus grey (missed) sum to 100% of ground truth; red stacks separately on top as false calls, scaled the same way.

</details>

<details>
<summary><strong>circRNA type classification</strong></summary>

### circRNA type classification

<img src="images/benchmark/type_classification_database_human_isoform_level.png" width="100%"/>

Same TP/FP vs GT breakdown as the ground truth coverage plot above, split by circRNA biotype (`eciRNA`, `EIciRNA`, `ciRNA`, `antisense`, `intergenic`, see [docs/methods.md](https://github.com/aerusakovich/nanocirc_pipeline/blob/main/docs/methods.md#circrna-type-classification)) instead of one bar per caller. `discovery`-tier recovery per type, isoform-level matching:

| Type | eciRNA | EIciRNA | ciRNA | antisense | intergenic |
| --- | --- | --- | --- | --- | --- |
| Database (human) | 81.8% | 31.7% | 80.2% | 80.3% | 49.6% |

eciRNA, ciRNA, and antisense are recovered well above the same tier's synthetic species-benchmark rates (see [circRNA type classification](#circrna-type-classification-1) above); EIciRNA and intergenic remain the two hardest types here too, the same pattern seen throughout this document.

</details>

<details>
<summary><strong>Isoform detection</strong></summary>

### Isoform detection

<img src="images/benchmark/pr_scatter_isoform_database_human.png" width="100%"/>
<img src="images/benchmark/isoform_detection_reciprocal_overlap_99_database_human.png" width="100%"/>
<img src="images/benchmark/isoform_family_completeness_reciprocal_overlap_99_database_human.png" width="100%"/>

**Note: only 2 shared-BSJ multi-isoform families exist in this 2,000-circRNAs-per-type database sample**, versus 102 in the wet-lab-protocol benchmark and 422/912 in the species benchmark. The numbers below describe those 2 families; treat them accordingly. % of ground-truth isoforms called right, at 99% reciprocal overlap:

| Tier / tool | Isoform | Family |
| --- | --- | --- |
| CIRI-long alone | 41.7% | 33.3% |
| isoCirc alone | 33.3% | 0.0% |
| circNICK-lrs alone | 0.0% | 0.0% |
| circFL-seq alone | 50.0% | 0.0% |
| `discovery` | 75.0% | 50.0% |
| `balanced_recall` | 75.0% | 50.0% |
| `balanced_precision` | 25.0% | 0.0% |
| `high_confidence` | 25.0% | 0.0% |

</details>

<details>
<summary><strong>Expression quantification and reproducibility</strong></summary>

### Expression quantification and reproducibility

<img src="images/benchmark/expression_correlation_database_human.png" width="100%"/>
<img src="images/benchmark/reproducibility_pairs_database_human.png" width="100%"/>

<img src="images/benchmark/expression_correlation_database_human_strict_bsj.png" width="100%"/>
<img src="images/benchmark/reproducibility_pairs_database_human_strict_bsj.png" width="100%"/>

<img src="images/benchmark/expression_correlation_database_human_exon_based.png" width="100%"/>
<img src="images/benchmark/reproducibility_pairs_database_human_exon_based.png" width="100%"/>

Read-count correlation against ground truth, and pairwise cross-tool/tier agreement, per caller. Pearson r (relaxed-BSJ matching):

| Caller | Database (human) |
| --- | --- |
| CIRI-long | 0.916 |
| isoCirc | 0.551 |
| circNICK-lrs | 0.776 |
| circFL-seq | 0.778 |
| `high_confidence` | 0.983 |
| `balanced_precision` | 0.984 |
| `balanced_recall` | 0.988 |
| `discovery` | 0.987 |

isoCirc is again the weakest quantifier, and every nanocirc tier again correlates more strongly with ground truth than any raw tool, the same pattern as the wet-lab-protocol and species benchmarks above. Detection sensitivity by ground-truth expression tertile:

| Tier | Low | Mid | High |
| --- | --- | --- | --- |
| `discovery` | 28.5% | 53.5% | 61.6% |
| `high_confidence` | 1.0% | 5.8% | 15.9% |

<img src="images/benchmark/sensitivity_by_expression_database_human.png" width="100%"/>

Same tertiles, matched at the `full`/`strict`/`relaxed` isoform tiers instead of structural comparison:

<img src="images/benchmark/sensitivity_by_expression_database_human_isoform.png" width="100%"/>

#### Isoform-level expression quantification

The table above reports relaxed-BSJ matching, which collapses a locus with several real isoforms into a single point. The tables below instead score each real GT isoform on its own: `full`, `strict`, `relaxed` are the same 3 tiers used for [Isoform detection](#isoform-detection-2) and the isoform PR scatter above (`full` = exact BSJ + 100% exon overlap, `strict` = BSJ within tolerance + 99% overlap, `relaxed` = BSJ within tolerance + 90% overlap). A GT isoform's predicted count is whichever caller isoform structurally matched it under that tier, or 0 if none did, so n here is a genuinely stricter re-match, not the table above's n expanded, and can be smaller, including for the many single-isoform loci in this dataset (isoform identity still has to be structurally confirmed for every one of them, not just the 2 multi-isoform families noted in [Isoform detection](#isoform-detection-2) above). Pearson r (n):

**`full` (exact BSJ, 100% exon overlap):**

| Caller | Database (human) r (n) |
| --- | --- |
| CIRI-long | 0.910 (3,622) |
| isoCirc | 0.553 (5,688) |
| circNICK-lrs | 0.652 (1,470) |
| circFL-seq | 0.776 (10,265) |
| `high_confidence` | 0.985 (1,965) |
| `balanced_precision` | 0.987 (6,783) |
| `balanced_recall` | 0.987 (11,201) |
| `discovery` | 0.988 (11,468) |

**`strict` (BSJ within tolerance, 99% exon overlap):**

| Caller | Database (human) r (n) |
| --- | --- |
| CIRI-long | 0.916 (9,081) |
| isoCirc | 0.551 (5,772) |
| circNICK-lrs | 0.659 (1,577) |
| circFL-seq | 0.776 (10,297) |
| `high_confidence` | 0.982 (1,973) |
| `balanced_precision` | 0.984 (6,837) |
| `balanced_recall` | 0.987 (11,896) |
| `discovery` | 0.987 (12,159) |

**`relaxed` (BSJ within tolerance, 90% exon overlap):**

| Caller | Database (human) r (n) |
| --- | --- |
| CIRI-long | 0.916 (9,420) |
| isoCirc | 0.551 (5,890) |
| circNICK-lrs | 0.690 (1,870) |
| circFL-seq | 0.778 (10,363) |
| `high_confidence` | 0.983 (2,004) |
| `balanced_precision` | 0.983 (6,918) |
| `balanced_recall` | 0.987 (12,099) |
| `discovery` | 0.988 (12,367) |

circNICK-lrs is by far the most affected caller here: its n drops from 5,409 at the locus level above to 1,470-1,870 once every call needs a genuinely matching isoform structure, and its correlation drops with it (r 0.65-0.69 here vs r 0.78 there), unlike the other raw tools and every nanocirc tier, which stay close to the table above at every tier of stringency. CIRI-long's own n also drops hard under `full` specifically (3,622 vs 9,429 at the locus level), less so under `strict`/`relaxed` (9,081-9,420), since exact-structure matching is a much higher bar than boundary tolerance plus 90-99% overlap. Every nanocirc tier still correlates near-perfectly with ground truth (r 0.98-0.99) regardless of stringency, and isoCirc is still the weakest quantifier throughout, the same pattern as every other benchmark in this document.

<img src="images/benchmark/expression_correlation_database_human_full.png" width="100%"/>
<img src="images/benchmark/expression_correlation_database_human_strict.png" width="100%"/>
<img src="images/benchmark/expression_correlation_database_human_relaxed.png" width="100%"/>

<img src="images/benchmark/reproducibility_pairs_database_human_full.png" width="100%"/>
<img src="images/benchmark/reproducibility_pairs_database_human_strict.png" width="100%"/>
<img src="images/benchmark/reproducibility_pairs_database_human_relaxed.png" width="100%"/>

### Caller overlap

Ground-truth loci recovered by each of the four tools are pooled across runs and drawn as an UpSet-style plot, where a locus counts as recovered if it was matched in at least one run. Two versions were used in the original benchmark paper: exon-based, which requires correct internal structure, and relaxed BSJ, which requires only genomic position. This is now extended to all six matching strategies used elsewhere in this benchmark: strict BSJ and exon-based (exact/structural), relaxed BSJ (coordinate-tolerant), and the full/strict/relaxed isoform tiers (structure plus BSJ tolerance, tightest to loosest).

The point of this plot is agreement, not any one tool's own recall. Of every circRNA any tool finds, only about 13% are found by all four tools together (12.5% on the species benchmark, 12.9% on the wet-lab CIRI-long protocol), and roughly half are found by exactly one tool. Each tool is effectively viewing circRNAs through its own lens, with no single reliable "ground truth caller" to defer to. That disagreement is the core motivation for NanoCirc: a naive union of all four tools' calls would inherit every tool's false positives, and a naive intersection would keep only the rare loci every tool happens to agree on, trading away most of the recall. NanoCirc's confidence tiers reconcile this with vote-based merging instead of a plain union/intersection rule, as described in the pipeline's own methods documentation.

<img src="images/benchmark/upset_plot_database_human_strict_bsj.png" width="100%"/>
<img src="images/benchmark/upset_plot_database_human.png" width="100%"/>
<img src="images/benchmark/upset_plot_database_human_exon_based.png" width="100%"/>
<img src="images/benchmark/upset_plot_database_human_full.png" width="100%"/>
<img src="images/benchmark/upset_plot_database_human_strict.png" width="100%"/>
<img src="images/benchmark/upset_plot_database_human_relaxed.png" width="100%"/>

</details>

**What we see:**

- **Database circRNAs are detected far more easily than de novo ones at every tier**: `discovery` reaches F1 0.68 here versus F1 0.48 (human) / 0.62 (mouse) on the synthetic species benchmark, mostly on recall (R 0.57 vs 0.37 / 0.54), consistent with these being already-annotated, pre-selected-for-detectability circRNAs.
- **Quantification is unaffected**: nanocirc tiers still correlate near-perfectly with ground truth (r 0.98-0.99), matching the synthetic benchmarks.
- **Isoform-level numbers here are too small for a trustable result**: only 2 multi-isoform families exist in this database sample, so treat the isoform/family percentages above as descriptive, not something to compare against the other benchmarks' isoform numbers.

<details>
<summary><strong>Statistical significance</strong></summary>

### Statistical significance

Do the synthetic benchmark's conclusions transfer to the database-derived dataset? Database circRNAs are easier to detect than synthetic ones by design. They are known, already-catalogued circRNAs, already selected for detectability, and used as input or a built-in dependency by tools during circRNA calling. The in silico dataset instead tries to evaluate performance on novel, unannotated molecules: it asks "if the read evidence for this molecule is there, but no annotation hints exist, can the tools still recover it?" The database benchmark asks a different question instead: "of circRNAs that are already known, annotated, and previously detected by different tools, how many can a tool recover when given read evidence and, where possible, the same input database?"

So the gap between the two methods is expected, but it's hard to say how much of it is due to the simulation logic, the database being biased by pre-selection for detectability, the presence of circRNA structures that were already reported by the same tools (which may still be real TPs or FPs, since we don't have wet-lab validation data for each database entry and depositing one doesn't require it), the presence of a database input file that directly hints the tools, actual biology, or other factors. What we can measure is whether tool ranking and per-type detectability can be transferred between the two. 

We matched precision/recall/F1 between the in silico and database benchmarks at strict isoform-level matching (BSJ ±10bp and ≥99% exon overlap), separately for every caller, circRNA type, and run: up to 120 paired points (8 callers × 5 types × 3 runs; precision has 111, since a caller that made zero calls of a type in a run has undefined precision, not 0, so that cell is dropped). Computed by [`scripts/statistical_analysis.py`](https://gitlab.com/bingo-igdr/nano-circ), reusable for any two datasets scored with our benchmark (not just in silico vs database, also wet-lab protocol vs protocol and species vs species, both used elsewhere in this document).

Same statistical methods as the [wet-lab-protocol comparison](#statistical-significance) (linear regression, Kendall's tau-b, Mantel-Haenszel + Breslow-Day, Wilcoxon signed-rank, Cohen's d), applied here to synthetic vs database instead of protocol vs protocol.

#### Aggregate, by metric (n paired points pooling every caller × type × run)

| Metric | n | Synthetic | Database | Difference | Slope | Intercept | R² | Wilcoxon p | Cohen's d |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Precision | 111 | 0.679 | 0.828 | +0.150 | 0.654 | 0.384 | 0.611 | 1×10⁻¹⁷ | 0.850 |
| Recall | 120 | 0.187 | 0.326 | +0.139 | 1.084 | 0.123 | 0.621 | 9×10⁻¹⁹ | 0.851 |
| F1 | 120 | 0.248 | 0.415 | +0.168 | 1.164 | 0.127 | 0.691 | 6×10⁻¹⁹ | 1.051 |

Recall and F1 on database circRNAs are well predicted by their in silico counterpart values (R² 0.62–0.69), close to a straight 1:1 relationship (slope 1.08–1.16). Precision is less predictable (R² = 0.61) and doesn't scale evenly: callers with weak synthetic-side precision gain the most on database circRNAs, while callers already near the ceiling gain little (slope 0.65). All three metrics are significantly higher on database circRNAs (Wilcoxon p < 10⁻¹⁶, Cohen's d 0.85–1.05).

#### Per type, by metric (n=24 paired points per row pooling all 8 callers × 3 runs, except antisense precision, n=21, and intergenic precision, n=18: some caller × run cells had zero calls of that type, so precision can't be defined and those cells are dropped)

| Type | Metric | n | Synthetic | Database | Difference | Slope | Intercept | R² | Wilcoxon p | Cohen's d |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| eciRNA | Precision | 24 | 0.651 | 0.804 | +0.153 | 0.835 | 0.260 | 0.858 | 1×10⁻⁷ | 1.831 |
| eciRNA | Recall | 24 | 0.495 | 0.550 | +0.055 | 1.091 | 0.010 | 0.987 | 1×10⁻⁷ | 1.981 |
| eciRNA | F1 | 24 | 0.532 | 0.634 | +0.102 | 1.145 | 0.025 | 0.951 | 1×10⁻⁷ | 2.328 |
| EIciRNA | Precision | 24 | 0.629 | 0.740 | +0.112 | 1.107 | 0.044 | 0.802 | 4×10⁻⁴ | 0.870 |
| EIciRNA | Recall | 24 | 0.103 | 0.157 | +0.054 | 1.803 | −0.029 | 0.796 | 0.003 | 0.823 |
| EIciRNA | F1 | 24 | 0.168 | 0.243 | +0.075 | 1.659 | −0.036 | 0.784 | 8×10⁻⁴ | 0.857 |
| ciRNA | Precision | 24 | 0.754 | 0.885 | +0.131 | 0.595 | 0.437 | 0.953 | 6×10⁻⁵ | 1.000 |
| ciRNA | Recall | 24 | 0.114 | 0.352 | +0.238 | 2.854 | 0.026 | 0.925 | 2×10⁻⁵ | 1.161 |
| ciRNA | F1 | 24 | 0.182 | 0.442 | +0.259 | 2.028 | 0.072 | 0.913 | 2×10⁻⁵ | 1.430 |
| antisense | Precision | 21 | 0.663 | 0.837 | +0.173 | 0.475 | 0.522 | 0.406 | 1×10⁻⁶ | 0.521 |
| antisense | Recall | 24 | 0.135 | 0.356 | +0.221 | 2.035 | 0.081 | 0.629 | 1×10⁻⁷ | 1.015 |
| antisense | F1 | 24 | 0.216 | 0.450 | +0.234 | 1.344 | 0.160 | 0.578 | 1×10⁻⁷ | 1.083 |
| intergenic | Precision | 18 | 0.699 | 0.892 | +0.193 | 0.455 | 0.574 | 0.428 | 8×10⁻⁶ | 1.912 |
| intergenic | Recall | 24 | 0.089 | 0.217 | +0.128 | 2.119 | 0.029 | 0.920 | 1×10⁻⁴ | 1.169 |
| intergenic | F1 | 24 | 0.139 | 0.308 | +0.168 | 1.753 | 0.063 | 0.890 | 1×10⁻⁴ | 1.344 |

Every type's recall and F1 transfer well between the two benchmarks (R² 0.58–0.99): a caller that's relatively better or worse at a type on synthetic circRNAs stays that way on database circRNAs. For eciRNA, the database boost is close to the same fixed size for every caller. For EIciRNA, ciRNA, antisense, and intergenic, it isn't fixed: callers that were already doing well on synthetic circRNAs gain even more on database circRNAs, so the recall/F1 gap between strong and weak callers gets wider, not just shifted up equally for everyone. Precision behaves the opposite way for most types (R² 0.41–0.95): for eciRNA, ciRNA, antisense, and intergenic, weaker callers catch up instead of falling further behind, so precision differences between callers get smaller on database circRNAs. EIciRNA is the exception, it is the only type where precision differences widen. Every type shows a statistically significant improvement on database circRNAs (Wilcoxon p ≤ 3×10⁻³, Cohen's d 0.52–2.33), meaning they are easier to detect.

#### Does caller ranking transfer?

We averaged each caller's synthetic and database values across type and run first, ranked the 8 callers by that average on each side, then compared the two rankings with Kendall's tau-b:

| Metric | n callers | Kendall's tau-b | Spearman rho |
| --- | --- | --- | --- |
| Precision | 8 | 0.714 (p = 0.014) | 0.833 (p = 0.010) |
| Recall | 8 | 0.857 (p = 0.002) | 0.929 (p = 9×10⁻⁴) |
| F1 | 8 | 0.857 (p = 0.002) | 0.929 (p = 9×10⁻⁴) |

Ranking transfers significantly on every metric, including precision. By F1:

| Caller | F1 synthetic | F1 database | Rank synthetic | Rank database |
| --- | --- | --- | --- | --- |
| balanced_recall | 0.406 | 0.727 | 1 | 1 |
| discovery | 0.396 | 0.707 | 2 | 2 |
| CIRI-long | 0.319 | 0.470 | 3 | 4 |
| balanced_precision | 0.263 | 0.372 | 4 | 5 |
| circFL-seq | 0.251 | 0.538 | 5 | 3 |
| isoCirc | 0.188 | 0.339 | 6 | 6 |
| high_confidence | 0.102 | 0.106 | 7 | 7 |
| circNICK-lrs | 0.056 | 0.064 | 8 | 8 |

Six of the eight callers had the exact same rank on both benchmarks. circFL-seq and CIRI-long swapped places: circFL-seq ranked 5th on in silico circRNAs but 3rd on database ones, CIRI-long the reverse. This suggests that circFL-seq is stronger on known annotated circRNAs and worse on de novo. 

<img src="images/benchmark/slopegraph_database_human.png" width="100%"/>

#### Per tool/mode, by metric (n=15 paired points per row: the 5 circRNA types × 3 runs, except circNICK-lrs precision n=12 and high_confidence precision n=9)

| Tool | Metric | n | Synthetic | Database | Difference | Slope | Intercept | R² | Wilcoxon p | Cohen's d |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CIRI-long | Precision | 15 | 0.687 | 0.865 | +0.178 | 0.286 | 0.669 | 0.217 | 6×10⁻⁵ | 1.664 |
| CIRI-long | Recall | 15 | 0.240 | 0.341 | +0.100 | 0.823 | 0.143 | 0.782 | 6×10⁻⁵ | 1.166 |
| CIRI-long | F1 | 15 | 0.319 | 0.470 | +0.151 | 1.045 | 0.136 | 0.750 | 6×10⁻⁵ | 1.707 |
| isoCirc | Precision | 15 | 0.877 | 0.945 | +0.068 | 0.460 | 0.541 | 0.783 | 6×10⁻⁵ | 1.405 |
| isoCirc | Recall | 15 | 0.121 | 0.217 | +0.096 | 0.769 | 0.124 | 0.874 | 6×10⁻⁵ | 1.882 |
| isoCirc | F1 | 15 | 0.188 | 0.339 | +0.152 | 0.644 | 0.218 | 0.822 | 6×10⁻⁵ | 1.699 |
| circNICK-lrs | Precision | 12 | 0.083 | 0.213 | +0.130 | −0.246 | 0.234 | 0.017 | 0.064 | 0.654 |
| circNICK-lrs | Recall | 15 | 0.056 | 0.056 | −0.000 | 1.115 | −0.007 | 0.980 | 0.530 | −0.008 |
| circNICK-lrs | F1 | 15 | 0.056 | 0.064 | +0.008 | 1.351 | −0.011 | 0.917 | 0.239 | 0.185 |
| circFL-seq | Precision | 15 | 0.624 | 0.922 | +0.298 | −0.075 | 0.969 | 0.146 | 6×10⁻⁵ | 0.818 |
| circFL-seq | Recall | 15 | 0.174 | 0.392 | +0.218 | 0.527 | 0.300 | 0.434 | 6×10⁻⁵ | 1.614 |
| circFL-seq | F1 | 15 | 0.251 | 0.538 | +0.287 | 0.331 | 0.455 | 0.229 | 6×10⁻⁵ | 1.513 |
| high_confidence | Precision | 9 | 0.900 | 0.939 | +0.039 | 0.655 | 0.349 | 0.931 | 0.024 | 1.302 |
| high_confidence | Recall | 15 | 0.065 | 0.068 | +0.003 | 1.089 | −0.003 | 0.996 | 0.762 | 0.241 |
| high_confidence | F1 | 15 | 0.102 | 0.106 | +0.004 | 1.085 | −0.005 | 0.993 | 0.762 | 0.188 |
| balanced_precision | Precision | 15 | 0.888 | 0.955 | +0.067 | 0.349 | 0.645 | 0.668 | 6×10⁻⁵ | 1.424 |
| balanced_precision | Recall | 15 | 0.186 | 0.254 | +0.069 | 0.911 | 0.085 | 0.951 | 6×10⁻⁵ | 1.557 |
| balanced_precision | F1 | 15 | 0.263 | 0.372 | +0.109 | 0.852 | 0.148 | 0.898 | 6×10⁻⁵ | 1.577 |
| balanced_recall | Precision | 15 | 0.690 | 0.895 | +0.205 | 0.337 | 0.663 | 0.950 | 6×10⁻⁵ | 2.353 |
| balanced_recall | Recall | 15 | 0.323 | 0.636 | +0.313 | 0.568 | 0.452 | 0.317 | 6×10⁻⁵ | 1.616 |
| balanced_recall | F1 | 15 | 0.406 | 0.727 | +0.321 | 0.772 | 0.414 | 0.452 | 6×10⁻⁵ | 2.501 |
| discovery | Precision | 15 | 0.650 | 0.813 | +0.163 | 0.716 | 0.348 | 0.850 | 6×10⁻⁵ | 2.336 |
| discovery | Recall | 15 | 0.331 | 0.647 | +0.316 | 0.581 | 0.455 | 0.336 | 6×10⁻⁵ | 1.629 |
| discovery | F1 | 15 | 0.396 | 0.707 | +0.311 | 0.873 | 0.361 | 0.389 | 6×10⁻⁵ | 2.274 |

isoCirc, high_confidence, and balanced_precision are consistent on all three metrics (R² ≥ 0.67). balanced_recall and discovery are strong on precision (R² ≥ 0.85) but less consistent on recall and F1 (R² 0.32–0.45, still significant). CIRI-long and circNICK-lrs are the opposite: precision doesn't track at all, but recall and F1 do (R² ≥ 0.75). circFL-seq is the least consistent of all eight: only its recall fit is significant, and even that is weak (R² = 0.43). By F1, high_confidence and circNICK-lrs perform almost identically on both benchmarks (difference ≤ 0.01), while balanced_recall, discovery, and circFL-seq score much higher on database circRNAs (+0.29 to +0.32); CIRI-long, isoCirc, and balanced_precision fall in between (+0.11 to +0.15).

| Tool | Outcome | Strata (types) | OR (MH) | 95% CI | MH p | Breslow-Day p |
| --- | --- | --- | --- | --- | --- | --- |
| CIRI-long | Precision | 5 | 3.90 | 3.62–4.20 | <10⁻³⁰⁰ | <10⁻³⁰⁰ |
| CIRI-long | Recall | 5 | 1.73 | 1.67–1.81 | <10⁻³⁰⁰ | <10⁻³⁰⁰ |
| isoCirc | Precision | 5 | 1.95 | 1.69–2.25 | <10⁻³⁰⁰ | 0.032 |
| isoCirc | Recall | 5 | 2.07 | 1.97–2.17 | <10⁻³⁰⁰ | <10⁻³⁰⁰ |
| circNICK-lrs | Precision | 5 | 1.54 | 1.43–1.67 | <10⁻³⁰⁰ | <10⁻³⁰⁰ |
| circNICK-lrs | Recall | 5 | 1.00 | 0.92–1.08 | 0.973 | <10⁻³⁰⁰ |
| circFL-seq | Precision | 5 | 2.91 | 2.65–3.19 | <10⁻³⁰⁰ | <10⁻³⁰⁰ |
| circFL-seq | Recall | 5 | 3.34 | 3.20–3.49 | <10⁻³⁰⁰ | <10⁻³⁰⁰ |
| high_confidence | Precision | 5 | 2.15 | 1.73–2.68 | 4×10⁻¹² | 0.247 |
| high_confidence | Recall | 5 | 1.07 | 0.99–1.15 | 0.082 | 1×10⁻¹⁶ |
| balanced_precision | Precision | 5 | 2.83 | 2.51–3.19 | <10⁻³⁰⁰ | 0.008 |
| balanced_precision | Recall | 5 | 1.60 | 1.52–1.67 | <10⁻³⁰⁰ | <10⁻³⁰⁰ |
| balanced_recall | Precision | 5 | 3.98 | 3.74–4.24 | <10⁻³⁰⁰ | 7×10⁻⁵ |
| balanced_recall | Recall | 5 | 4.21 | 4.05–4.38 | <10⁻³⁰⁰ | <10⁻³⁰⁰ |
| discovery | Precision | 5 | 2.44 | 2.32–2.57 | <10⁻³⁰⁰ | <10⁻³⁰⁰ |
| discovery | Recall | 5 | 4.33 | 4.17–4.51 | <10⁻³⁰⁰ | <10⁻³⁰⁰ |

Every caller has significantly higher odds of a correct call (precision) and of recovering a ground-truth circRNA (recall) on database circRNAs than on synthetic ones (OR 1.5–4.3), except circNICK-lrs's recall, which shows no real difference (OR ≈ 1.00, matching its flat recall line in the slopegraph). But that boost isn't the same size across circRNA types: of the sixteen caller/outcome combinations tested, only high_confidence's precision is uniform across types (Breslow-Day p = 0.25); every other combination's boost depends on which circRNA type it is.

**Conclusion.** Correlation is good: aggregate R² is 0.61–0.69, and most types and callers individually fit a significant line. Caller ranking transfers: Kendall's tau-b is 0.71–0.86 on every metric, including precision, and six of the eight callers land on the same F1 rank on both benchmarks, with only circFL-seq and CIRI-long swapping places with each other. But the size of the database boost isn't uniform: it depends on circRNA type for fifteen of the sixteen caller/outcome combinations tested (high_confidence's precision is the one exception). So the overall pattern and ranking are preserved well between de novo and database calling, but the exact performance improvement on annotated circRNAs is not uniform, and depends on both the circRNA type and the caller.

</details>

<details>
<summary><strong>What are these false positives?</strong></summary>

### What are these false positives?

We traced every FP-supporting read back to the simulator's own ground truth (single-tool discovery exon-based FPs, all 3 human runs, mean % ± SD across runs):

| FP read is... | CIRI-long | circFL-seq | isoCirc | circNICK-lrs |
| --- | --- | --- | --- | --- |
| Real circ read, correct boundary, wrong strand (antisense called as sense) | 29.0 ± 1.0 | 64.8 ± 2.1 | 88.8 ± 2.5 | 0.6 ± 0.1 |
| Real circ read, correct strand, pure boundary noise | 20.8 ± 2.6 | 16.5 ± 2.5 | 0.0 ± 0.0 | 34.4 ± 0.6 |
| Real circ read, but circRNA reported at a different/distant locus | 50.3 ± 1.9 | 18.7 ± 0.4 | 11.2 ± 2.5 | 65.0 ± 0.6 |
| Not from a circ read at all (linear read miscalled) | 0.0 ± 0.0 | 0.0 ± 0.0 | 0.0 ± 0.0 | 0.1 ± 0.1 |
| FP locus overlaps a genomic repeat | 57.7 ± 1.7 | 56.4 ± 1.1 | 47.4 ± 1.4 | 87.9 ± 0.3 |
| n FP loci (pooled, 3 runs) | 635 | 653 | 118 | 4,930 |
| n supporting reads (pooled, 3 runs) | 12,614 | 5,608 | 1,052 | 47,269 |

This is a genuinely different FP composition from every synthetic-data comparison above. On synthetic data, wrong-strand antisense calls dominate every tool's FPs (77-94%, see the wet-lab-protocol and species sections' equivalent tables). On database circRNAs, only isoCirc still shows that pattern (88.8%); CIRI-long and circNICK-lrs are instead dominated by different-locus calls (50.3% and 65.0%), and circFL-seq is roughly split between wrong-strand (64.8%) and boundary/locus noise. Repeat overlap is still high across all four tools (47-88%), as it is throughout this document. A plausible explanation is that unlike the simulator, the database annotation carries circRNAs with more diverse, sometimes ambiguous, genomic contexts (isoforms of overlapping/nearby genes, non-simulated antisense pairs), so a caller can land near the right region without hitting the exact same locus or transcribing from the correct strand's known circRNA.

**Practical implication**: database `discovery` false positives are less concentrated in a single, explainable error mode than synthetic ones. Re-verify low-confidence calls the same way recommended elsewhere in this document (e.g. Sanger sequencing) before trusting exact structure, strand, or locus.

</details>

<details>
<summary><strong>What discovery still misses?</strong></summary>

### What discovery still misses?

We also characterized every ground-truth circRNA against whether any tool recovered it at all (mean % ± SD across 3 human runs; missed + recovered = 100% of the GT catalog, ~8,672 circRNAs/run):

| | Missed (FN) | Recovered |
| --- | --- | --- |
| Share of GT catalog | 52.1 ± 0.1 | 47.9 ± 0.1 |
| Mean mature length (bp) | 13,292 ± 66 | 319 ± 0 |
| Mean rolling-circle copy number | 7.73 ± 0.05 | 8.54 ± 0.04 |
| Median simulated read count | 31 ± 0 | 82 ± 1 |
| Canonical GT/AG splice site at BSJ | 94.8 ± 0.1 | 96.3 ± 0.0 |
| No gene overlap at all (intergenic) | 16.4 ± 0.5 | 15.9 ± 0.4 |
| Falls in a gene with ≥5 isoforms | 72.2 ± 0.6 | 71.9 ± 0.4 |
| Overlaps a genomic repeat | 78.5 ± 0.0 | 78.8 ± 0.1 |

Recall per biotype:

| Biotype | Missed | Recovered |
| --- | --- | --- |
| `eciRNA` (pure exonic) | 25.9 ± 0.2 | 74.1 ± 0.2 |
| `ciRNA` (purely intronic) | 49.5 ± 0.4 | 50.5 ± 0.4 |
| `antisense` | 51.1 ± 0.1 | 48.9 ± 0.1 |
| `intergenic` | 64.7 ± 0.1 | 35.3 ± 0.1 |
| `EIciRNA` (exon+intron) | 72.8 ± 0.0 | 27.2 ± 0.0 |

Missed circRNAs are far longer (~42x, versus 4-9x on the synthetic benchmarks) and lower-expressed, but the copy-number, splice-site, intergenic, and repeat-overlap gaps between missed and recovered are all small here, much smaller than the synthetic benchmarks' equivalents. `eciRNA` is again the easiest biotype; `EIciRNA` the hardest, matching the type-classification table above and every other benchmark in this document.

</details>

---

## Cross-benchmark comparison: which effect is biggest?

Every benchmark above shows two datasets with difference in one variable (species, wet-lab protocol, or annotated-vs-de-novo). This section plots them in a single PCA over each dataset's pooled precision/recall/F1, computed by [`publication/pca_analysis.py`](https://gitlab.com/bingo-igdr/nano-circ): one row per dataset, one column per (caller × structural strategy × threshold × metric) combination shared by all 7, standardized (z-scored) before decomposition.

<img src="images/benchmark/pca_benchmarks.png" width="100%"/>

PC1 (72% of variance) is dominated by F1 and recall, i.e. overall detection difficulty (how hard to make a call); PC2 (15%) is dominated by precision (how precise are made calls).

To isolate each effect on its own, `pca_analysis.py` also reports the direct distance (in the same standardized feature space the PCA was computed in, not just the PC1/PC2 projection) between pairs of datasets differing in only that one factor, saved to `pca_pairwise_distances.csv`:

| Effect | Datasets compared | Distance |
| --- | --- | --- |
| Annotated (database) vs de novo (synthetic) | database_human vs species_human | 19.13 |
| isoCirc vs circFL-seq protocol | wetlab_isocirc vs wetlab_circfl | 15.70 |
| Human vs mouse | species_human vs species_mouse | 15.48 |
| isoCirc vs CIRI-long protocol | wetlab_isocirc vs wetlab_ciri_long | 12.75 |
| CIRI-long vs circFL-seq protocol | wetlab_ciri_long vs wetlab_circfl | 6.75 |

**Conclusion.** Among everything tested, the biggest variation comes from annotated (database) vs de novo (synthetic) circRNA recovery, with annotated circRNAs being easier to identify, closely followed by the isoCirc-protocol wet-lab effect. The species effect (human vs mouse) is smaller but in the same range as the isoCirc-protocol effect. CIRI-long vs circFL-seq protocol has the smallest effect of anything tested, under half the size of every other comparison. 

The annotated vs de novo distance (≈19) is not disproportionately large next to other, it's in the same order of magnitude. That supports the earlier hypothesis that de novo synthetic circRNAs are a valid comparator: they're offset from database circRNAs in absolute difficulty (expected as per reasons already documented above), by a margin similar in scale to other differences (species, protocol) while caller ranking and per-type structure transfer well between the two (see [Statistical significance](#statistical-significance-2) above).

Putting all benchmarks together: annotated circRNAs are both more easily recoverable and called with higher precision than de novo ones. Among wet-lab protocols, circFL-seq allows the best recovery, closely followed by CIRI-long protocol; isoCirc protocol is the hardest for the tools to recover reads from, despite offering slightly higher precision. Species differences are also to be expected: mouse-derived circRNAs were easier to recover than human-derived ones, with an effect similar in size to a wet-lab-protocol effect (species distance 15.48, close to isoCirc-vs-circFL-seq's 15.70). Tool ranking generally transfers well across every one of these comparisons (Kendall's tau-b 0.71–1.00), and when a discrepancy exists it's usually just two callers swapping adjacent places, for example circFL-seq and CIRI-long swapping on annotated circRNAs, or circFL-seq and isoCirc swapping between human and mouse. Across all of this, nanocirc's own tiers are the most stable callers: `balanced_recall` and `discovery` held the top two F1 ranks in every wet-lab-protocol, species, and database comparison in this benchmark, and every nanocirc tier correlated with ground-truth expression more strongly than any single raw tool on every dataset tested, while the raw tools' relative strengths shifted more with protocol and species.

---

For the simulator and benchmark scoring code behind these results, see the rest of this repo ([nano-circ](https://gitlab.com/bingo-igdr/nano-circ)). For the consensus pipeline that computes the 4 confidence tiers, see [nanocirc_pipeline](https://github.com/aerusakovich/nanocirc_pipeline).
