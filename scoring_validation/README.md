# SIREN: Scoring System Validation

This report documents how the the SIREN scoring weights were validated by grid search using the script validate_siren_designVIII_FULL_v4.py, and how to interpret the outputs in designVIII_FULL_report_topk20/. Scoring weights were stress-tested using an existing benchmark folder of completed SIREN runs, and what the resulting grid search says about **redundancy control**, **region stability**, and **which weight settings are effectively equivalent**.

All results below correspond to running:

```bash
python validate_siren_designVIII_FULL_v4.py \
  --bench_dir=bench_results \
  --out_dir=designVIII_FULL_report_topk20 \
  --topk=20 \
  --rnai_lengths=100,150,200,250,300,350,400 \
  --step=4 \
  --p_single='-0.005,-0.01,-0.02,-0.03,-0.05,-0.07,-0.1,-0.15,-0.2' \
  --p_extra='0,-0.1,-0.5,-1,-2,-3,-5,-7,-10,-15,-20,-30,-40,-50,-60' \
  --baseline_ps=-0.1 --baseline_pe=-30 \
  --tol=0.001 \
  --ratio_bins=30 \
  --make_plots
```

The core machine-readable outputs referenced here are `grid.tsv`, `best.tsv`, and `runs.tsv` written inside `designVIII_FULL_report_topk20/`

---

## What was validated

The goal of this validation is to test whether the scoring weights used by **designVIII** behave as intended across many real SIREN outputs:

1. **Redundancy control:** penalizing multiple distinct siRNAs that hit the same off-target (the “redundancy penalty”) should reduce redundant off-targeting in the *selected* dsRNA regions.
2. **Region stability:** changing scoring weights should *not* cause erratic changes in which regions are selected, except in regimes where the scoring is qualitatively different (e.g., removing redundancy penalty entirely).
3. **Parameter robustness / plateau:** if many different weight pairs produce essentially the same results, that’s a useful property (it means the system does not depend on a fragile “magic number”).

This validation is performed on **real SIREN run outputs** in `bench_results/` rather than on synthetic toy data.

---

## Benchmark data used

`bench_results/` is treated as a collection of independent SIREN runs (one “run folder” per target). For each run folder, the validator expects the same core products that SIREN generates, especially:

- the target sequence (`target.fa`), and  
- the summarized off-target hits (`off_targets_summary.tsv`).

From `runs.tsv` in your report archive, the validation aggregated results across **464 run folders**, with per-length run counts:

- L=100: 464 runs  
- L=150: 464 runs  
- L=200: 464 runs  
- L=250: 459 runs  
- L=300: 455 runs  
- L=350: 450 runs  
- L=400: 446 runs  

(Those small drops at larger lengths are expected in practice: some targets are too short or otherwise fail length-specific requirements.)

---

## How candidate dsRNA regions were generated

For each run folder and each requested dsRNA length **L**, the script enumerates candidate dsRNA windows along the target sequence using a fixed stride:

- `--step=4` means the candidate window **start coordinate advances by 4 nt** each time.
- This yields a consistent, reproducible candidate set per target and per length.

For each candidate region, the validator computes scoring features from the run’s off-target summary: essentially, what off-target interactions are implied by the siRNAs contained within that dsRNA window.

---

## The weight grid that was evaluated

Each scoring system is defined by two penalties:

- `p_single`: penalty applied when an off-target is hit (first time).
- `p_extra`: additional penalty applied for **redundancy**, i.e., extra siRNAs that also hit the *same* off-target beyond the first one.

The script evaluates the full Cartesian product:

- 9 values of `p_single`
- 15 values of `p_extra`

That is **135 scoring systems per length**, applied across hundreds of runs.

A specific pair is designated as the reference baseline:

- `baseline_ps=-0.1`
- `baseline_pe=-30`

---

## What the main metrics mean

These are computed over the **Top-K selected dsRNA regions** (here `--topk=20`), then averaged across runs.

- **`mean_unique_topk`**: the average number of **distinct off-target entities** engaged by the Top-K selected dsRNA windows. This behaves like a “breadth” measure: higher values mean the selected regions still produce many unique off-target interactions (not redundancy per se, just the number of unique off-targets hit at least once).

- **`mean_extra_topk`**: a redundancy-focused measure quantifying **how much additional off-targeting occurs because multiple siRNAs within a region hit the same off-target**. This is the main “redundancy burden” signal and the primary objective minimized in `best.tsv`.

- **`mean_total_topk`**: the overall off-targeting load implied by the Top-K regions (magnitude of off-target interactions), without restricting to “unique only”.

Separately, the validator computes region-agreement metrics **against the baseline** (baseline is used only as a reference point for overlap, not for redundancy magnitude):

- **`jaccard_starts_vs_baseline`**: Jaccard overlap of the **set of start positions** of Top-K windows vs the baseline Top-K.
- **`jaccard_cov_vs_baseline`**: Jaccard overlap of **nucleotide coverage** of Top-K windows vs baseline Top-K.
- The fraction of runs whose Top-K changes vs baseline is summarized as `frac_runs_topk_changed`.

The remaining CLI knobs:

- `--tol=0.001` is used only for plateau detection: values within **0.1%** of the best objective are treated as effectively equivalent.
- `--ratio_bins=30` bins the penalty ratio into 30 bins for ratio-trend plots.

---

## Key result 1: the “best” weights are identical to the baseline (up to scaling)

From `grid.tsv` and `best.tsv`, the best-selected configuration for every dsRNA length was:

- `p_single = -0.2`, `p_extra = -60` (ratio = 300)

and this produces the same aggregated metrics as the baseline:

- baseline: `p_single = -0.1`, `p_extra = -30` (ratio = 300)

Because designVIII scoring is a linear sum of penalties, multiplying both penalties by the same constant preserves ranking. That’s why (−0.1, −30) and (−0.2, −60) yield identical Top-K choices and identical summary metrics across all lengths in this grid.

This directly supports the interpretation that, for this benchmark, the penalty **ratio** is the meaningful knob, and the absolute magnitudes mainly rescale scores without changing rankings.

---

## Key result 2: removing the redundancy penalty produces a dramatic outlier regime

A direct ablation confirms the redundancy term is doing real work:

- ablation: `p_extra = 0`, with `p_single = -0.1`

When redundancy is removed, redundancy metrics explode and the selected regions change massively relative to the baseline. A compact snapshot (Top-K = 20, step = 4; values are averages across runs):

| dsRNA length | baseline mean_extra_topk | ablation mean_extra_topk | baseline mean_total_topk | ablation mean_total_topk | ablation Jaccard (starts vs baseline) | fraction of runs with Top-K start set changed |
|---:|---:|---:|---:|---:|---:|---:|
| 100 | 471.97 | 2741.18 | 983.50 | 4255.14 | 0.217 | 0.963 |
| 200 | 1464.86 | 7102.64 | 2919.42 | 9673.22 | 0.269 | 0.972 |
| 300 | 2908.99 | 12856.26 | 5471.43 | 16506.63 | 0.275 | 0.965 |
| 400 | 5409.77 | 21737.11 | 9236.29 | 26212.68 | 0.274 | 0.956 |

Two things matter here.

First, `mean_extra_topk` increases by roughly **4×–6×** when redundancy is not penalized, and `mean_total_topk` rises sharply too.

Second, the region selection becomes fundamentally different: roughly **95–97%** of runs select a different Top-K start set compared to the baseline and overlap collapses (start-set Jaccard ~0.22–0.28). This is strong internal evidence that the redundancy term is not cosmetic; it directly prevents selection of regions that generate many redundant off-targeting siRNAs.

---

## Key result 3: once redundancy penalty is non-zero, performance quickly plateaus

Outside the `p_extra=0` regime, the grid shows a broad plateau where many non-zero weight pairs behave similarly in aggregate. This is why untransformed heatmaps can look “flat”: one outlier regime dominates the dynamic range, and the remaining structure is subtle.

The log-compressed heatmaps in the output folder are intended to make that subtle structure visible rather than being dominated by the `p_extra=0` outlier.

---

## Key result 4: redundancy increases strongly with dsRNA length (even under baseline)

Under the baseline scoring (or any equivalent ratio-300 scaling), the absolute redundancy burden grows rapidly with dsRNA length:

- L=100: `mean_extra_topk ≈ 472`
- L=200: `mean_extra_topk ≈ 1465`
- L=300: `mean_extra_topk ≈ 2909`
- L=400: `mean_extra_topk ≈ 5410`

This is expected: longer dsRNAs generate more candidate siRNAs, increasing the opportunity for repeated hits on the same off-target.

---

## Practical conclusion for defaults

Within the tested grid and this benchmark dataset:

- The default weights **(−0.1, −30)** are already on the optimum plateau for redundancy minimization in this experiment.
- The grid’s “best” solution (−0.2, −60) is just a scaled version of the baseline and therefore does not justify changing defaults.
- Removing redundancy penalty (`p_extra=0`) creates a clearly inferior regime with dramatically higher redundancy and unstable region selection.

A simple recommendation consistent with these results is: keep **(−0.1, −30)** as the default, and interpret the **ratio** as the meaningful hyperparameter; absolute scaling mainly rescales scores.

---

## Where to look in the output folder

To reproduce or independently verify:

- `grid.tsv`: aggregated metrics for every `(p_single, p_extra, rnai_len)` combination.
- `best.tsv`: best configuration per length (and overall) under the selected objective.
- `runs.tsv`: per-run/per-config summaries that underlie the aggregate statistics.
- Plots:
  - redundancy heatmaps (`hm_log1p_redundancy_len*.png`)
  - region overlap heatmaps (`hm_jaccard_starts_len*.png`, `hm_jaccard_cov_len*.png`)
  - ratio summaries (`ratio_curve_*.png`, `heatmap_min_redundancy_by_ratio.png`)

