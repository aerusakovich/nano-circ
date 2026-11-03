from __future__ import annotations

import argparse
import logging
import sys

from .config import CIRC_TYPES, SimConfig


def _add_simulate_parser(subparsers):
    p = subparsers.add_parser("simulate", help="Simulate circRNAs from a genome + GTF")
    p.add_argument("--genome", required=True, help="Reference genome FASTA (indexed or indexable)")
    p.add_argument("--gtf", required=True, help="Gene annotation GTF")
    p.add_argument("--outdir", required=True, help="Output directory")
    p.add_argument("--config", help="YAML config (defaults ship in configs/default.yaml if omitted)")
    p.add_argument("--total", type=int, help="Total circRNAs to simulate, apportioned via config.type_fractions")
    for t in CIRC_TYPES:
        p.add_argument(f"--n-{t}", type=int, default=None, help=f"Exact number of {t} circRNAs (overrides --total apportionment)")
    p.add_argument("--total-reads", type=int, default=None, help="Sequencing depth budget for n_reads_target allocation")
    p.add_argument("--isoform-mean", type=float, default=None, help="Average distinct isoforms drawn per used locus")
    p.add_argument("--isoform-max", type=int, default=None, help="Cap on isoforms per locus")
    p.add_argument("--min-reads-per-circrna", type=int, default=None,
                    help="Every circRNA gets at least this many reads targeted, so low-TPM circRNAs aren't left with zero")
    p.add_argument("--min-circle-length", type=int, default=None,
                    help="Reject circRNAs shorter than this (too short to ever produce a valid long read)")
    p.add_argument("--max-circle-length", type=int, default=None,
                    help="Reject circRNAs longer than this (no realistic long read could span the junction)")
    p.add_argument("--seed", type=int, default=None, help="Random seed (overrides config.seed)")
    p.add_argument("--quiet", action="store_true", help="Suppress per-type progress bars")
    return p


def _add_fit_config_parser(subparsers):
    p = subparsers.add_parser("fit-config", help="Fit a config YAML from a real circRNA database BED12")
    p.add_argument("--db-bed", required=True, help="Real circRNA database BED12 (e.g. circbase/circatlas merged)")
    p.add_argument("--genome", required=True, help="Reference genome FASTA matching --db-bed coordinates")
    p.add_argument("--gtf", required=True, help="Gene annotation GTF matching --db-bed coordinates")
    p.add_argument("--out", required=True, help="Path to write the fitted config YAML")
    p.add_argument("--base-config", help="Config to use as a base for fields dbfit can't estimate (defaults included)")
    return p


def _add_check_bsj_parser(subparsers):
    p = subparsers.add_parser("check-bsj", help="Mandatory BSJ-in-read verification for simulated/real long reads")
    p.add_argument("--fastq", required=True, help="FASTQ of long reads to check")
    p.add_argument("--mature-fasta", required=True, help="FASTA of mature (non-concatemer) circRNA sequences")
    p.add_argument("--metadata", required=True, help="circRNA_metadata.tsv or ground_truth.tsv from `simulate`")
    p.add_argument("--outdir", required=True, help="Output directory")
    p.add_argument("--minimap2", default="minimap2", help="Path to minimap2 executable")
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--write-filtered-fastq", action="store_true",
                    help="Also write filtered_confirmed_only.fastq: linear reads kept as-is, "
                         "circular-labeled reads kept only if they actually crossed the BSJ")
    return p


def _add_topup_plan_parser(subparsers):
    p = subparsers.add_parser(
        "topup-plan", help="List circRNAs with too little BSJ-confirmed read support, for a follow-up read-sim round")
    p.add_argument("--ground-truth", required=True, help="ground_truth_reconciled.tsv from `check-bsj`")
    p.add_argument("--min-confirmed", type=int, default=1, help="Flag circRNAs with fewer confirmed reads than this")
    p.add_argument("--out", required=True, help="Path to write the NanoSim-compatible top-up expression profile")
    return p


def _add_pool_select_parser(subparsers):
    p = subparsers.add_parser(
        "pool-select",
        help="From a one-shot over-generated mega-pool, select up to n_reads_target BSJ-confirmed reads per circRNA")
    p.add_argument("--bsj-reads", required=True, help="bsj_check_reads.tsv from `check-bsj` on the mega-pool")
    p.add_argument("--metadata", required=True, help="circRNA_metadata.tsv (source of each circRNA's n_reads_target)")
    p.add_argument("--mega-fastq", required=True, help="The over-generated mega-pool FASTQ that was BSJ-checked")
    p.add_argument("--outdir", required=True, help="Output directory")
    p.add_argument("--subsample-to", type=int, default=None,
                    help="Draw one replicate run of approximately this many reads (Poisson per circRNA)")
    p.add_argument("--seed", type=int, default=0,
                    help="Seed for --subsample-to's per-circRNA Poisson draw and read shuffle")
    p.add_argument("--ground-truth", default=None,
                    help="ground_truth.tsv from `simulate`. With --subsample-to, also writes "
                         "run_circrnas.bed and run_ground_truth.tsv for this run")
    p.add_argument("--shortfall-profile-out", default=None,
                    help="Also write a NanoSim expression profile (TPM = target - available) "
                         "for every circRNA short of its target, for a follow-up top-up round")
    return p


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="nanocirc-sim", description="Fast, reproducible, database-fittable circRNA simulator")
    parser.add_argument("-v", "--verbose", action="store_true", help="Debug logging")
    subparsers = parser.add_subparsers(dest="command", required=True)
    _add_simulate_parser(subparsers)
    _add_fit_config_parser(subparsers)
    _add_check_bsj_parser(subparsers)
    _add_topup_plan_parser(subparsers)
    _add_pool_select_parser(subparsers)
    return parser


def _load_config(args) -> SimConfig:
    cfg = SimConfig.from_yaml(args.config) if args.config else SimConfig.default()
    if args.seed is not None:
        cfg.seed = args.seed
    if args.total_reads is not None:
        cfg.total_reads = args.total_reads
    if args.isoform_mean is not None:
        cfg.isoform_mean = args.isoform_mean
    if args.isoform_max is not None:
        cfg.isoform_max = args.isoform_max
    if args.min_reads_per_circrna is not None:
        cfg.min_reads_per_circrna = args.min_reads_per_circrna
    if args.min_circle_length is not None:
        cfg.min_circle_length = args.min_circle_length
    if args.max_circle_length is not None:
        cfg.max_circle_length = args.max_circle_length
    return cfg


def _cmd_simulate(args) -> int:
    from .simulate import run_simulation

    cfg = _load_config(args)
    explicit_counts = {t: getattr(args, f"n_{t}") for t in CIRC_TYPES}
    result = run_simulation(
        genome_path=args.genome, gtf_path=args.gtf, config=cfg, outdir=args.outdir,
        total=args.total, explicit_counts=explicit_counts, seed=cfg.seed,
        show_progress=not args.quiet,
    )
    stats = result["type_stats"]
    print("\nPer-type results:")
    for t in CIRC_TYPES:
        if t in stats:
            s = stats[t]
            print(f"  {t:12s} requested={s['requested']:6d} achieved={s['achieved']:6d} "
                  f"attempts={s['attempts']:6d} rejection={s['rejection_rate']*100:.1f}%")
    print(f"\nTotal circRNAs: {len(result['candidates'])}")
    print(f"Elapsed: {result['elapsed_seconds']:.1f}s")
    print(f"Outputs written to: {args.outdir}")
    return 0


def _cmd_fit_config(args) -> int:
    from .dbfit import fit_config_from_database

    base_cfg = SimConfig.from_yaml(args.base_config) if args.base_config else SimConfig.default()
    fitted = fit_config_from_database(
        db_bed_path=args.db_bed, genome_path=args.genome, gtf_path=args.gtf, base_config=base_cfg,
    )
    fitted.to_yaml(args.out)
    print(f"Fitted config written to: {args.out}")
    return 0


def _cmd_check_bsj(args) -> int:
    from .bsj_check import check_bsj

    summary = check_bsj(
        fastq_path=args.fastq, mature_fasta_path=args.mature_fasta, metadata_path=args.metadata,
        outdir=args.outdir, minimap2_path=args.minimap2, threads=args.threads,
        write_filtered=args.write_filtered_fastq,
    )
    print(f"BSJ check complete: {summary['n_confirmed']}/{summary['n_circ_labeled']} "
          f"circular-labeled reads confirmed to actually cross the BSJ "
          f"({summary['n_confirmed'] / max(summary['n_circ_labeled'], 1) * 100:.1f}%)")
    print(f"Outputs written to: {args.outdir}")
    return 0


def _cmd_topup_plan(args) -> int:
    from .topup import select_deficient_circles, write_topup_expression_profile

    deficient = select_deficient_circles(args.ground_truth, min_confirmed=args.min_confirmed)
    write_topup_expression_profile(deficient, args.out)
    print(f"{len(deficient)} circRNA(s) with fewer than {args.min_confirmed} confirmed BSJ reads")
    print(f"Top-up expression profile written to: {args.out}")
    return 0


def _cmd_pool_select(args) -> int:
    import os

    from .pool_select import (
        deficient_circle_ids, filter_plan_to_selected, plan_run_targets, plan_selection,
        run_selected_read_ids, select_run_reads, selected_read_ids, undersupplied_shortfalls,
        write_run_bed, write_run_ground_truth_tsv, write_selected_fastq, write_selection_stats,
        write_shortfall_profile,
    )

    os.makedirs(args.outdir, exist_ok=True)
    plan = plan_selection(args.bsj_reads, args.metadata)

    # Computed from the full plan, not any one --subsample-to run.
    deficient = deficient_circle_ids(plan)
    with open(os.path.join(args.outdir, "deficient_circle_ids.txt"), "w") as fh:
        fh.write("\n".join(deficient) + ("\n" if deficient else ""))

    if args.shortfall_profile_out:
        shortfalls = undersupplied_shortfalls(plan)
        write_shortfall_profile(shortfalls, args.shortfall_profile_out)
        print(f"{len(shortfalls)} circRNA(s) short of their target "
              f"(total shortfall {sum(shortfalls.values())} reads), profile written to: "
              f"{args.shortfall_profile_out}")

    if args.subsample_to is not None:
        run_targets = plan_run_targets(plan, args.subsample_to, args.seed)
        run_reads = select_run_reads(plan, run_targets, args.seed)
        selected = run_selected_read_ids(run_reads)
        stats_plan = filter_plan_to_selected(plan, selected)
        if args.ground_truth:
            n_bed = write_run_bed(args.ground_truth, run_targets, os.path.join(args.outdir, "run_circrnas.bed"))
            print(f"Wrote {n_bed} circRNAs to run_circrnas.bed (score = reads picked for this run)")
            n_tsv = write_run_ground_truth_tsv(
                args.ground_truth, run_targets, os.path.join(args.outdir, "run_ground_truth.tsv"))
            print(f"Wrote {n_tsv} circRNAs to run_ground_truth.tsv (full columns, same per-run filtering)")
    else:
        selected = selected_read_ids(plan)
        stats_plan = plan

    stats = write_selection_stats(stats_plan, os.path.join(args.outdir, "selection_stats.tsv"))
    n_written = write_selected_fastq(args.mega_fastq, selected, os.path.join(args.outdir, "selected_reads.fastq"))

    print(f"Selected {n_written} reads across {stats['n_circrnas']} circRNAs "
          f"({stats['n_deficient']} with zero confirmed reads available in the mega-pool)")
    print(f"Outputs written to: {args.outdir}")
    return 0


def main(argv=None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO,
                         format="%(asctime)s - %(levelname)s - %(message)s")

    handlers = {
        "simulate": _cmd_simulate, "fit-config": _cmd_fit_config, "check-bsj": _cmd_check_bsj,
        "topup-plan": _cmd_topup_plan, "pool-select": _cmd_pool_select,
    }
    return handlers[args.command](args)


if __name__ == "__main__":
    sys.exit(main())
