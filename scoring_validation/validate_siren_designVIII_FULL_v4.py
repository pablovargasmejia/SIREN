#!/usr/bin/env python3
"""
validate_siren_designVIII_FULL.py

ONE script that does EVERYTHING, end-to-end:
  - Reads bench_results run folders
  - Recomputes siren_designVIII-like window scoring from off_targets_summary.tsv + target.fa
  - Sweeps p_single / p_extra and rnai_length values
  - Computes redundancy + total hits + unique hits in TopK
  - Computes region-similarity between scoring systems (TopK start Jaccard + base-coverage Jaccard)
  - Picks BEST scoring system (per length + overall across lengths)
  - Computes "ratio plateau" (where |p_extra|/|p_single| saturates)
  - Writes ALL TSVs + ALL plots into a single output folder
  - Writes a short Markdown report you can paste into rebuttal

Inputs required per run folder (inside --bench_dir):
  - off_targets_summary.tsv
  - other_files/target.fa   (or target.fa)

Key definitions (matches your design logic)
-------------------------------------------
For each candidate dsRNA window, define for each off-target transcript t:
  count_t = number of siRNAs from transcript t fully contained in the dsRNA window

Redundancy (extra hits) = sum_t max(0, count_t - 1)

We compute it exactly using:
  total_hits  = sum_t count_t
  targets_hit = number of t with count_t > 0
  extra_hits  = total_hits - targets_hit     (exact identity)

Scoring function (same structure as designVIII):
  score = p_single * unique_hits + p_extra * extra_hits
(higher score is better; penalties are negative)

"Best" objective you requested
------------------------------
"best regions = less redundancy and less siRNAs (unique or redundant)"

So we rank parameter sets by:
  1) minimize mean_extra_topk   (redundancy)
  2) then minimize mean_total_topk
  3) then minimize mean_unique_topk
  4) then maximize mean_jaccard_cov_vs_baseline (tie-break for stability)

Plateau (ratio saturation)
--------------------------
For each rnai_len and metric M:
  M* = minimum over the full grid
  best_M_at_or_above_ratio(r) = min{ M : ratio_abs_pe_ps >= r }
Saturation ratio = smallest r such that best_M_at_or_above_ratio(r) <= (1+tol)*M*

Usage (copy/paste)
------------------
cd /fs/scratch/PAS1755/Pablo/SIREN/siren17

python validate_siren_designVIII_FULL.py \
  --bench_dir bench_results \
  --out_dir designVIII_FULL_report \
  --topk 10 \
  --rnai_lengths 200,250,300 \
  --step 4 \
  --p_single '-0.005,-0.01,-0.02,-0.03,-0.05,-0.07,-0.1,-0.15,-0.2' \
  --p_extra  '0,-0.1,-0.5,-1,-2,-3,-5,-7,-10,-15,-20,-30,-40,-50,-60' \
  --baseline_ps -0.1 --baseline_pe -30 \
  --tol 0.001 \
  --ratio_bins 30 \
  --make_plots

Outputs are all inside --out_dir:
  grid.tsv, runs.tsv, best.tsv, plateau_thresholds.tsv, ratio_curves_*.tsv
  many PNG plots
  report.md
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import List, Optional, Tuple

import numpy as np
import pandas as pd

def df_to_markdown_table(df: pd.DataFrame) -> str:
    """Minimal markdown table renderer (no external deps)."""
    if df is None or len(df) == 0:
        return "_(empty)_"
    sdf = df.copy()
    # stringify + keep readable floats
    for c in sdf.columns:
        def _fmt(x):
            if pd.isna(x):
                return ""
            if isinstance(x, float):
                return f"{x:.6g}"
            return str(x)
        sdf[c] = sdf[c].map(_fmt)
    cols = list(sdf.columns)
    widths = {c: max(len(c), *(sdf[c].map(len).tolist())) for c in cols}

    def fmt_row(vals):
        return "| " + " | ".join(v.ljust(widths[c]) for v, c in zip(vals, cols)) + " |"

    header = fmt_row(cols)
    sep = "| " + " | ".join("-" * widths[c] for c in cols) + " |"
    body = "\n".join(fmt_row([row[c] for c in cols]) for _, row in sdf.iterrows())
    return "\n".join([header, sep, body])



SIRNA_COORD_RE = re.compile(r"(\d+)-(\d+)")


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bench_dir", default="bench_results")
    ap.add_argument("--out_dir", default="designVIII_FULL_report")
    ap.add_argument("--topk", type=int, default=10)
    ap.add_argument("--rnai_lengths", default="200,300")
    ap.add_argument("--step", type=int, default=4)

    ap.add_argument("--p_single", required=True)
    ap.add_argument("--p_extra", required=True)

    ap.add_argument("--baseline_ps", type=float, default=-0.1)
    ap.add_argument("--baseline_pe", type=float, default=-30.0)

    ap.add_argument("--tol", type=float, default=0.001)
    ap.add_argument("--ratio_bins", type=int, default=30)

    ap.add_argument("--max_runs", type=int, default=0)
    ap.add_argument("--make_plots", action="store_true")
    return ap.parse_args()


def parse_list_floats(s: str) -> List[float]:
    return [float(x.strip()) for x in str(s).split(",") if x.strip()]


def parse_list_ints(s: str) -> List[int]:
    return [int(x.strip()) for x in str(s).split(",") if x.strip()]


def read_fasta_len(fp: Path) -> int:
    n = 0
    with fp.open("r") as f:
        for line in f:
            if line.startswith(">"):
                continue
            n += len(line.strip())
    return n


def find_target_fa(run_dir: Path) -> Optional[Path]:
    for c in [
        run_dir / "other_files" / "target.fa",
        run_dir / "target.fa",
        run_dir / "other_files" / "target.fasta",
        run_dir / "target.fasta",
    ]:
        if c.exists():
            return c
    return None


def sirna_coords(name: str) -> Tuple[int, int]:
    m = list(SIRNA_COORD_RE.finditer(str(name)))
    if not m:
        raise ValueError(f"Cannot parse siRNA coords from: {name}")
    m = m[-1]
    return int(m.group(1)), int(m.group(2))


def load_offtargets(run_dir: Path) -> Optional[pd.DataFrame]:
    fp = run_dir / "off_targets_summary.tsv"
    if not fp.exists():
        return None
    df = pd.read_csv(fp, sep="\t")
    if "siRNA names" not in df.columns:
        raise ValueError(f"{fp} missing column 'siRNA names'")
    return df


def build_arrays(off_targets: pd.DataFrame):
    """
    Returns arrays for fast per-window metrics.
    """
    pairs = []
    all_starts = []
    per_target_unique_starts = []

    for names_str in off_targets["siRNA names"].astype(str).tolist():
        items = [x.strip() for x in names_str.split(",") if x.strip()]
        starts_set = set()
        for it in items:
            s, e = sirna_coords(it)
            pairs.append((s, e))
            all_starts.append(s)
            starts_set.add(s)
        per_target_unique_starts.append(np.array(sorted(starts_set), dtype=np.int64))

    if pairs:
        lens = [e - s + 1 for (s, e) in pairs]
        sirna_len = int(pd.Series(lens).mode().iloc[0])
        unique_starts = np.array(sorted({s for (s, e) in pairs}), dtype=np.int64)
    else:
        sirna_len = 21
        unique_starts = np.array([], dtype=np.int64)

    all_starts = np.array(sorted(all_starts), dtype=np.int64)
    return unique_starts, all_starts, per_target_unique_starts, sirna_len


def idx_interval_from_start_range(a: int, b: int, step: int, max_start: int) -> Optional[Tuple[int, int]]:
    """
    Candidate starts: 1, 1+step, 1+2*step, ...
    Map coordinate interval [a,b] into candidate index interval [i0,i1].
    """
    if b < 1 or a > max_start:
        return None
    a = max(a, 1)
    b = min(b, max_start)
    if a > b:
        return None
    i0 = (a - 1 + (step - 1)) // step  # ceil
    i1 = (b - 1) // step              # floor
    if i0 > i1:
        return None
    return int(i0), int(i1)


def intervals_for_target(starts_unique: np.ndarray, W: int, max_start: int, step: int):
    """
    A window start s contains a siRNA start p if p in [s, s+W-1] (containment start bound).
    So s in [p-W+1, p]. Convert to index intervals and merge.
    """
    if starts_unique.size == 0:
        return []
    idx_intervals = []
    for p in starts_unique:
        a = int(p) - W + 1
        b = int(p)
        iv = idx_interval_from_start_range(a, b, step=step, max_start=max_start)
        if iv is not None:
            idx_intervals.append(iv)
    if not idx_intervals:
        return []
    idx_intervals.sort()
    merged = [idx_intervals[0]]
    for a, b in idx_intervals[1:]:
        la, lb = merged[-1]
        if a <= lb + 1:
            merged[-1] = (la, max(lb, b))
        else:
            merged.append((a, b))
    return merged


def compute_features_for_length(target_len: int, rnai_len: int, step: int,
                                unique_starts: np.ndarray, all_starts: np.ndarray,
                                per_target_unique_starts, sirna_len: int):
    """
    Candidate windows are spaced by --step, and siRNAs must be fully contained.
    W = rnai_len - sirna_len + 1 is the number of valid siRNA start positions in a window.
    """
    max_start = target_len - rnai_len + 1
    if max_start <= 0:
        return None
    W = rnai_len - sirna_len + 1
    if W <= 0:
        return None

    starts = np.arange(1, max_start + 1, step, dtype=np.int64)
    end_allowed = starts + W - 1  # == win_end - (sirna_len-1)

    # unique_hits
    if unique_starts.size == 0:
        unique_hits = np.zeros_like(starts, dtype=np.int64)
    else:
        lu = np.searchsorted(unique_starts, starts, side="left")
        ru = np.searchsorted(unique_starts, end_allowed, side="right")
        unique_hits = (ru - lu).astype(np.int64)

    # total_hits
    if all_starts.size == 0:
        total_hits = np.zeros_like(starts, dtype=np.int64)
    else:
        la = np.searchsorted(all_starts, starts, side="left")
        ra = np.searchsorted(all_starts, end_allowed, side="right")
        total_hits = (ra - la).astype(np.int64)

    # targets_hit
    n_cand = starts.size
    diff = np.zeros(n_cand + 1, dtype=np.int64)
    for t_starts in per_target_unique_starts:
        merged = intervals_for_target(t_starts, W=W, max_start=max_start, step=step)
        for i0, i1 in merged:
            diff[i0] += 1
            if i1 + 1 < diff.size:
                diff[i1 + 1] -= 1
    targets_hit = np.cumsum(diff)[:n_cand]

    # redundancy
    extra_hits = np.maximum(total_hits - targets_hit, 0)

    return starts, unique_hits, total_hits, targets_hit, extra_hits


def topk_idx(scores: np.ndarray, k: int) -> np.ndarray:
    kk = min(k, scores.size)
    return np.argsort(scores)[::-1][:kk]


def jaccard_set(a: np.ndarray, b: np.ndarray) -> float:
    A = set(map(int, a.tolist()))
    B = set(map(int, b.tolist()))
    if not A and not B:
        return 1.0
    if not A or not B:
        return 0.0
    return len(A & B) / len(A | B)


def merge_intervals(intervals: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
    if not intervals:
        return []
    intervals = sorted(intervals)
    merged = [intervals[0]]
    for s, e in intervals[1:]:
        ls, le = merged[-1]
        if s <= le + 1:
            merged[-1] = (ls, max(le, e))
        else:
            merged.append((s, e))
    return merged


def jaccard_coverage(intsA: List[Tuple[int, int]], intsB: List[Tuple[int, int]]) -> float:
    A = merge_intervals(intsA)
    B = merge_intervals(intsB)
    if not A and not B:
        return 1.0
    if not A or not B:
        return 0.0
    i = j = 0
    inter = 0
    while i < len(A) and j < len(B):
        a1, a2 = A[i]
        b1, b2 = B[j]
        s = max(a1, b1)
        e = min(a2, b2)
        if s <= e:
            inter += (e - s + 1)
        if a2 < b2:
            i += 1
        else:
            j += 1
    def total_len(M):
        return sum(e - s + 1 for s, e in M)
    union = total_len(A) + total_len(B) - inter
    return inter / union if union > 0 else 1.0


def best_at_or_above_ratio(sub: pd.DataFrame, ratio_thresholds: np.ndarray, metric: str) -> pd.DataFrame:
    out = []
    for r in ratio_thresholds:
        s2 = sub[sub["ratio_abs_pe_ps"] >= r]
        out.append((r, float(s2[metric].min()) if len(s2) else np.nan))
    return pd.DataFrame(out, columns=["ratio_threshold", f"best_{metric}_at_or_above_ratio"])


def make_plots(agg: pd.DataFrame, out_dir: Path):
    import matplotlib.pyplot as plt

    # Heatmaps per rnai_len on p_single x p_extra
    for L in sorted(agg["rnai_len"].unique()):
        sub = agg[agg["rnai_len"] == L].copy()

        def hm(value_col: str, title: str, outname: str, log1p: bool = False):
            piv = sub.pivot(index="p_single", columns="p_extra", values=value_col).sort_index().sort_index(axis=1)
            Z = piv.values
            if log1p:
                Z = np.log1p(Z)
            fig = plt.figure(figsize=(11, 6))
            ax = fig.add_subplot(111)
            im = ax.imshow(Z, aspect="auto", origin="lower")
            ax.set_title(title)
            ax.set_xlabel("p_extra")
            ax.set_ylabel("p_single")
            ax.set_xticks(np.arange(piv.shape[1]))
            ax.set_xticklabels([str(x) for x in piv.columns], rotation=90)
            ax.set_yticks(np.arange(piv.shape[0]))
            ax.set_yticklabels([str(y) for y in piv.index])
            fig.colorbar(im, ax=ax, label=(f"log1p({value_col})" if log1p else value_col))
            fig.tight_layout()
            fig.savefig(out_dir / outname, dpi=200)
            plt.close(fig)

        hm("mean_extra_topk", f"Total redundancy in TopK (mean_extra_topk) — rnai_len={L}", f"hm_total_redundancy_len{L}.png", log1p=False)
        hm("mean_extra_topk", f"log1p redundancy in TopK — rnai_len={L}", f"hm_log1p_redundancy_len{L}.png", log1p=True)
        hm("mean_total_topk", f"log1p total siRNA hits in TopK — rnai_len={L}", f"hm_log1p_totalhits_len{L}.png", log1p=True)
        hm("mean_unique_topk", f"log1p unique siRNA hits in TopK — rnai_len={L}", f"hm_log1p_uniquehits_len{L}.png", log1p=True)
        hm("mean_jaccard_cov_vs_baseline", f"Coverage Jaccard vs baseline — rnai_len={L}", f"hm_covJ_vs_baseline_len{L}.png", log1p=False)
        hm("mean_jaccard_starts_vs_baseline", f"Start Jaccard vs baseline — rnai_len={L}", f"hm_startJ_vs_baseline_len{L}.png", log1p=False)

    # Scatter: redundancy vs total hits
    fig = plt.figure(figsize=(9, 6))
    ax = fig.add_subplot(111)
    for L in sorted(agg["rnai_len"].unique()):
        sub = agg[agg["rnai_len"] == L]
        ax.scatter(sub["mean_total_topk"].to_numpy(), sub["mean_extra_topk"].to_numpy(), s=14, label=f"rnai_len={L}")
    ax.set_title("Tradeoff: total siRNA hits vs redundancy (TopK=10)")
    ax.set_xlabel("mean_total_topk")
    ax.set_ylabel("mean_extra_topk")
    ax.legend()
    fig.tight_layout()
    fig.savefig(out_dir / "scatter_total_vs_redundancy.png", dpi=200)
    plt.close(fig)

    # Scatter: ratio vs redundancy (log-x)
    fig = plt.figure(figsize=(9, 6))
    ax = fig.add_subplot(111)
    for L in sorted(agg["rnai_len"].unique()):
        sub = agg[agg["rnai_len"] == L]
        ax.scatter(sub["ratio_abs_pe_ps"].to_numpy(), sub["mean_extra_topk"].to_numpy(), s=14, label=f"rnai_len={L}")
    ax.set_xscale("log")
    ax.set_title("Grid points: redundancy vs ratio_abs_pe_ps")
    ax.set_xlabel("ratio_abs_pe_ps (log scale)")
    ax.set_ylabel("mean_extra_topk")
    ax.legend()
    fig.tight_layout()
    fig.savefig(out_dir / "scatter_ratio_vs_redundancy.png", dpi=200)
    plt.close(fig)


def make_ratio_analysis(agg: pd.DataFrame, out_dir: Path, tol: float, ratio_bins: int):
    import matplotlib.pyplot as plt

    metrics = ["mean_extra_topk", "mean_total_topk", "mean_unique_topk"]

    # saturation thresholds + curves
    rows = []
    for L in sorted(agg["rnai_len"].unique()):
        sub = agg[agg["rnai_len"] == L].dropna(subset=["ratio_abs_pe_ps"] + metrics).copy()
        ratios = np.array(sorted(sub["ratio_abs_pe_ps"].unique()), dtype=float)
        rmin = max(1e-9, float(np.min(ratios)))
        rmax = float(np.max(ratios))
        ratio_thresholds = np.unique(np.concatenate([ratios, np.logspace(np.log10(rmin), np.log10(rmax), 200)]))
        ratio_thresholds.sort()

        curve_df = pd.DataFrame({"ratio_threshold": ratio_thresholds})
        for m in metrics:
            m_star = float(sub[m].min())
            c = best_at_or_above_ratio(sub, ratio_thresholds, m)
            curve_df = curve_df.merge(c, on="ratio_threshold", how="left")
            ok = curve_df[f"best_{m}_at_or_above_ratio"] <= (1.0 + tol) * m_star
            r_sat = float(curve_df.loc[ok, "ratio_threshold"].iloc[0]) if ok.any() else float("nan")
            rows.append({"rnai_len": L, "metric": m, "min_value": m_star, "tol": tol, "ratio_saturation": r_sat, "max_ratio_tested": rmax})

        curve_df.to_csv(out_dir / f"ratio_curves_len{L}.tsv", sep="\t", index=False)

        # plot curve
        fig = plt.figure(figsize=(10, 6))
        ax = fig.add_subplot(111)
        for m in metrics:
            ax.plot(curve_df["ratio_threshold"].to_numpy(),
                    curve_df[f"best_{m}_at_or_above_ratio"].to_numpy(),
                    label=m)
        ax.set_xscale("log")
        ax.set_title(f"Best-achievable metrics vs ratio threshold (rnai_len={L})")
        ax.set_xlabel("ratio threshold (require ratio_abs_pe_ps >= r) [log]")
        ax.set_ylabel("best achievable value (lower is better)")
        ax.legend()
        fig.tight_layout()
        fig.savefig(out_dir / f"ratio_curves_len{L}.png", dpi=200)
        plt.close(fig)

    thr = pd.DataFrame(rows)
    thr.to_csv(out_dir / "ratio_plateau_thresholds.tsv", sep="\t", index=False)

    # Heatmap: rnai_len × ratio bins of min redundancy
    all_rat = agg["ratio_abs_pe_ps"].replace([np.inf, -np.inf], np.nan).dropna().to_numpy()
    rmin = max(1e-9, float(np.min(all_rat)))
    rmax = float(np.max(all_rat))
    edges = np.logspace(np.log10(rmin), np.log10(rmax), ratio_bins + 1)
    tmp = agg.copy()
    tmp["ratio_bin"] = pd.cut(tmp["ratio_abs_pe_ps"], bins=edges, include_lowest=True)

    hm = (
        tmp.groupby(["rnai_len", "ratio_bin"], observed=True)
           .agg(
               min_redundancy=("mean_extra_topk", "min"),
               min_total=("mean_total_topk", "min"),
               min_unique=("mean_unique_topk", "min"),
           )
           .reset_index()
    )
    piv = hm.pivot(index="rnai_len", columns="ratio_bin", values="min_redundancy")
    fig = plt.figure(figsize=(12, 4))
    ax = fig.add_subplot(111)
    im = ax.imshow(np.log1p(piv.values), aspect="auto", origin="lower")
    ax.set_title("log1p(min redundancy) by rnai_len × ratio bin")
    ax.set_xlabel("ratio_abs_pe_ps bins (log-spaced)")
    ax.set_ylabel("rnai_len")
    ax.set_yticks(np.arange(piv.shape[0]))
    ax.set_yticklabels([str(x) for x in piv.index])
    ax.set_xticks(np.arange(piv.shape[1]))
    ax.set_xticklabels([str(c) for c in piv.columns], rotation=90, fontsize=7)
    fig.colorbar(im, ax=ax, label="log(1+min mean_extra_topk)")
    fig.tight_layout()
    fig.savefig(out_dir / "heatmap_min_redundancy_by_ratio.png", dpi=200)
    plt.close(fig)

    return thr


def write_report(out_dir: Path, best_per_len: pd.DataFrame, best_overall: pd.DataFrame,
                 baseline: pd.DataFrame, ablation: pd.DataFrame, plateau: pd.DataFrame, topk: int, tol: float):
    md = []
    md.append("# SIREN designVIII scoring validation report\n")
    md.append(f"- TopK = **{topk}** windows per run\n")
    md.append(f"- Plateau tolerance = **{tol}** (relative)\n")

    md.append("\n## Best scoring parameters\n")
    md.append("### Best per RNAi length\n")
    md.append(df_to_markdown_table(best_per_len))
    md.append("\n\n### Best overall (regardless RNAi length)\n")
    md.append(df_to_markdown_table(best_overall))

    if len(baseline):
        md.append("\n\n## Baseline (-0.1, -30)\n")
        md.append(df_to_markdown_table(baseline))
    if len(ablation):
        md.append("\n\n## Ablation (p_extra=0, p_single=-0.1)\n")
        md.append(df_to_markdown_table(ablation))

    md.append("\n\n## Ratio plateau (saturation)\n")
    md.append("Smallest ratio |p_extra|/|p_single| where further increasing the ratio does not improve the metric beyond tolerance.\n")
    md.append(df_to_markdown_table(plateau))

    (out_dir / "report.md").write_text("\n".join(md))


def main():
    args = parse_args()
    bench_dir = Path(args.bench_dir).resolve()
    out_dir = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    rnai_lengths = parse_list_ints(args.rnai_lengths)
    ps_list = parse_list_floats(args.p_single)
    pe_list = parse_list_floats(args.p_extra)

    run_dirs = sorted([p for p in bench_dir.iterdir() if p.is_dir()])
    if args.max_runs and args.max_runs > 0:
        run_dirs = run_dirs[:args.max_runs]

    per_run_rows = []
    processed = 0
    for rd in run_dirs:
        ot = load_offtargets(rd)
        tf = find_target_fa(rd)
        if ot is None or tf is None:
            continue

        target_len = read_fasta_len(tf)
        unique_starts, all_starts, per_target_unique_starts, sirna_len = build_arrays(ot)

        for L in rnai_lengths:
            feats = compute_features_for_length(target_len, L, args.step,
                                                unique_starts, all_starts,
                                                per_target_unique_starts, sirna_len)
            if feats is None:
                continue

            starts, unique_hits, total_hits, targets_hit, extra_hits = feats

            # baseline topK for similarity
            base_scores = args.baseline_ps * unique_hits + args.baseline_pe * extra_hits
            base_idx = topk_idx(base_scores, args.topk)
            base_starts = starts[base_idx]
            base_ints = [(int(s), int(s + L - 1)) for s in base_starts]

            for ps in ps_list:
                for pe in pe_list:
                    scores = ps * unique_hits + pe * extra_hits
                    idx = topk_idx(scores, args.topk)
                    sel_starts = starts[idx]
                    sel_ints = [(int(s), int(s + L - 1)) for s in sel_starts]

                    per_run_rows.append({
                        "run_folder": rd.name,
                        "rnai_len": L,
                        "sirna_len_est": sirna_len,
                        "step": args.step,
                        "p_single": ps,
                        "p_extra": pe,
                        "ratio_abs_pe_ps": abs(pe) / abs(ps) if ps != 0 else float("inf"),
                        "mean_extra_topk": float(extra_hits[idx].mean()) if idx.size else float("nan"),
                        "mean_total_topk": float(total_hits[idx].mean()) if idx.size else float("nan"),
                        "mean_unique_topk": float(unique_hits[idx].mean()) if idx.size else float("nan"),
                        "jaccard_starts_vs_baseline": jaccard_set(sel_starts, base_starts),
                        "jaccard_cov_vs_baseline": jaccard_coverage(sel_ints, base_ints),
                    })

        processed += 1
        if processed % 50 == 0:
            print(f"Processed {processed} run folders...")

    if not per_run_rows:
        raise SystemExit("No usable runs found (need off_targets_summary.tsv + target.fa).")

    runs_df = pd.DataFrame(per_run_rows)

    # Aggregate across runs
    agg = runs_df.groupby(["rnai_len", "p_single", "p_extra", "ratio_abs_pe_ps"], as_index=False).agg(
        mean_extra_topk=("mean_extra_topk", "mean"),
        mean_total_topk=("mean_total_topk", "mean"),
        mean_unique_topk=("mean_unique_topk", "mean"),
        mean_jaccard_starts_vs_baseline=("jaccard_starts_vs_baseline", "mean"),
        mean_jaccard_cov_vs_baseline=("jaccard_cov_vs_baseline", "mean"),
        n_runs=("run_folder", "nunique"),
    )

    runs_df.to_csv(out_dir / "runs.tsv", sep="\t", index=False)
    agg.to_csv(out_dir / "grid.tsv", sep="\t", index=False)

    # BEST per rnai_len and overall (your requested objective)
    cols = ["rnai_len", "p_single", "p_extra", "ratio_abs_pe_ps",
            "mean_extra_topk", "mean_total_topk", "mean_unique_topk",
            "mean_jaccard_cov_vs_baseline", "mean_jaccard_starts_vs_baseline", "n_runs"]

    best_per_len = []
    for L in sorted(agg["rnai_len"].unique()):
        sub = agg[agg["rnai_len"] == L].copy()
        best = sub.sort_values(
            ["mean_extra_topk", "mean_total_topk", "mean_unique_topk", "mean_jaccard_cov_vs_baseline"],
            ascending=[True, True, True, False],
        ).head(1)
        best_per_len.append(best)
    best_per_len = pd.concat(best_per_len, ignore_index=True)

    best_overall = agg.sort_values(
        ["mean_extra_topk", "mean_total_topk", "mean_unique_topk", "mean_jaccard_cov_vs_baseline"],
        ascending=[True, True, True, False],
    ).head(1)

    best_tbl = pd.concat([best_per_len, best_overall.assign(note="BEST_OVERALL")], ignore_index=True)
    best_tbl.to_csv(out_dir / "best.tsv", sep="\t", index=False)

    baseline = agg[(agg.p_single == args.baseline_ps) & (agg.p_extra == args.baseline_pe)].sort_values("rnai_len")
    ablation = agg[(agg.p_single == args.baseline_ps) & (agg.p_extra == 0.0)].sort_values("rnai_len")

    # Ratio plateau
    plateau = make_ratio_analysis(agg, out_dir, tol=args.tol, ratio_bins=args.ratio_bins)

    # Plots
    if args.make_plots:
        make_plots(agg, out_dir)

    # Report
    write_report(out_dir, best_per_len[cols], best_overall[cols], baseline[cols], ablation[cols], plateau, args.topk, args.tol)

    # Print key tables
    print("\n=== BEST per rnai_len ===")
    print(best_per_len[cols].to_string(index=False))
    print("\n=== BEST overall (regardless rnai_len) ===")
    print(best_overall[cols].to_string(index=False))

    print("\nWrote everything into:")
    print(" ", out_dir)
    print("Key files:")
    print(" ", out_dir / "grid.tsv")
    print(" ", out_dir / "best.tsv")
    print(" ", out_dir / "ratio_plateau_thresholds.tsv")
    print(" ", out_dir / "report.md")


if __name__ == "__main__":
    main()
