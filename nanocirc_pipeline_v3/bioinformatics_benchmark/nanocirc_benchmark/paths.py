"""Dataset path resolution for one species within a benchmark dataset_root.

Directory layout expected under `dataset_root` (produced upstream by the
nf-core/nanocirc pipeline plus the 4 competing tools, external to this
package):

  fastq/<run>.fastq.gz
  ground_truth/<run>_ground_truth.tsv
  results/<species>/circrna/<run>/...                  per-run tool + NanoCirc outputs
  results/<species>/circrna/crossrun/<species>/<tier>/<species>_<tier>_crossrun.bed12
  benchmark/gt_prefiltered/<run>_gt_prefiltered.bed12   GT restricted to this run's expressed loci
  benchmark/<species>.genome                            bedtools-style chrom-sizes file
  benchmark/...                                         everything else is written by this package
"""
from pathlib import Path


class Paths:
    def __init__(self, dataset_root, species: str, runs: list, gtf_path=None):
        self.dataset_root = Path(dataset_root)
        self.species = species
        self.runs = runs
        # Only required by type_classification.py's circFL-style classifier.
        self.gtf_path = Path(gtf_path) if gtf_path is not None else None

    @property
    def classification_data_dir(self) -> Path:
        # Must be species-specific. gene.bed/exon.bed are built once and
        # reused (build_feature_beds skips rebuilding if already present).
        # Sharing this dir across species would silently reuse one species'
        # gene/exon annotation for another, and misclassify everything as
        # intergenic (no real coordinate overlap between the two genomes).
        d = self.bench_dir / "data" / "classification" / self.species
        d.mkdir(parents=True, exist_ok=True)
        return d

    # ---- top-level ----
    @property
    def fastq_dir(self) -> Path:
        return self.dataset_root / "fastq"

    @property
    def results_dir(self) -> Path:
        return self.dataset_root / "results" / self.species / "circrna"

    @property
    def bench_dir(self) -> Path:
        d = self.dataset_root / "benchmark"
        d.mkdir(parents=True, exist_ok=True)
        return d

    @property
    def plot_dir(self) -> Path:
        d = self.bench_dir / "plot"
        d.mkdir(parents=True, exist_ok=True)
        return d

    @property
    def bedtools_work_dir(self) -> Path:
        d = self.bench_dir / "bedtools_work"
        d.mkdir(parents=True, exist_ok=True)
        return d

    @property
    def genome_sizes_file(self) -> Path:
        path = self.bench_dir / f"{self.species}.genome"
        if not path.exists():
            raise FileNotFoundError(
                f"genome sizes file not found: {path}\n"
                "Expected a pre-built bedtools-style chrom-sizes file (chrom<TAB>length per line)."
            )
        return path

    # ---- per-run ----
    def run_fastq(self, run: str) -> Path:
        return self.fastq_dir / f"{run}.fastq.gz"

    def run_dir(self, run: str) -> Path:
        return self.results_dir / run

    def ground_truth_tsv(self, run: str) -> Path:
        return self.dataset_root / "ground_truth" / f"{run}_ground_truth.tsv"

    def gt_prefiltered_bed12(self, run: str) -> Path:
        return self.bench_dir / "gt_prefiltered" / f"{run}_gt_prefiltered.bed12"

    def nanocirc_clean_counts_tsv(self, run: str) -> Path:
        return self.run_dir(run) / "merged" / "clean" / f"{run}_discovery_clean_with_counts.tsv"

    # ---- crossrun-merged tier membership (group-level, not per-sample) ----
    def crossrun_tier_bed12(self, tier: str) -> Path:
        return self.results_dir / "crossrun" / self.species / tier / f"{self.species}_{tier}_crossrun.bed12"

    # ---- native per-tool raw output (per-read IDs / native BSJ_reads) ----
    def native_tool_paths(self, run: str) -> dict:
        base = self.run_dir(run)
        return {
            "isocirc": base / "isocirc" / "isocirc_output" / "isocirc.out",
            "circfl": base / "circfl_seq" / f"{run}_circfl" / "RG" / "circFL_Normal_pass.txt",
            "cirilong": base / "ciri_long" / f"{run}_cirilong" / f"{run}.reads",
            "circnick": base / "circnick_lrs" / f"{run}_circnick" / run / f"{run}.circRNA_candidates.annotated.txt",
        }

    # ---- structurally-merged BED12 per caller (tools + NanoCirc tiers),
    # the one convention every matching stage uses ----
    def caller_bed_path(self, run: str, key: str, tiers: list) -> Path:
        base = self.run_dir(run)
        if key == "cirilong":
            return base / "bed12" / f"{run}_cirilong.bed12"
        if key == "isocirc":
            return base / "isocirc" / "isocirc_output" / "isocirc.bed"
        if key == "circnick":
            return base / "bed12" / f"{run}_circnick.bed12"
        if key == "circfl":
            return base / "circfl_seq" / f"{run}_circfl" / "circFL_final.bed"
        if key in tiers:
            return base / "merged" / f"{run}_{key}.bed12"
        raise ValueError(f"unknown caller key: {key}")

    def caller_paths(self, run: str, tools: list, tiers: list) -> dict:
        return {key: self.caller_bed_path(run, key, tiers) for key in list(tools) + list(tiers)}
