#!/usr/bin/env python3
import argparse
import re
from pathlib import Path
from collections import Counter

import pandas as pd
from Bio import SeqIO


def norm(seq: str) -> str:
    return re.sub(r"\s+", "", str(seq).upper()).replace("U", "T")


def rc(seq: str) -> str:
    comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return norm(seq).translate(comp)[::-1]


def overlap(a_start: int, a_end: int, b_start: int, b_end: int) -> int:
    return max(0, min(a_end, b_end) - max(a_start, b_start))


def is_low_complexity(kmer: str, max_major_frac: float = 0.80) -> bool:
    kmer = kmer.upper()
    if not kmer:
        return True
    c = Counter(kmer)
    major = c.most_common(1)[0][1]
    if major / len(kmer) >= max_major_frac:
        return True
    # crude alternating-repeat rejection
    if len(set(kmer[::2])) <= 2 and len(set(kmer[1::2])) <= 2:
        return True
    return False


def find_design_tsv(run_dir: Path) -> Path | None:
    preferred = [
        run_dir / "rna_sequences_with_scores_and_primers.tsv",
        run_dir / "rna_sequences_with_scores.tsv",
        run_dir / "rna_sequences_with_scores_and_primers.txt",
        run_dir / "rna_sequences_with_scores.txt",
    ]
    for p in preferred:
        if p.exists() and p.stat().st_size > 0:
            return p

    for p in sorted(run_dir.glob("*.tsv")):
        name = p.name.lower()
        if ("rna" in name and "score" in name) or ("primer" in name and "rna" in name):
            if p.stat().st_size > 0:
                return p
    return None


def load_true_target_sequence(run_dir: Path) -> str | None:
    """
    CRITICAL: the real gene target is in other_files/target.fa (SIREN output structure).
    Never use targets_prefiltered.fa for target mapping (that is off-target DB).
    """
    candidates = [
        run_dir / "other_files" / "target.fa",
        run_dir / "target.fa",  # fallback
    ]
    for p in candidates:
        if p.exists() and p.stat().st_size > 0:
            rec = next(SeqIO.parse(str(p), "fasta"), None)
            if rec:
                return norm(str(rec.seq))
    return None


def detect_seq_col(df: pd.DataFrame) -> str:
    preferred = {
        "rnai_sequence", "rnai sequence", "sequence",
        "dsrna_sequence", "dsrna sequence", "rnai", "dsrna"
    }
    for c in df.columns:
        if c.lower() in preferred:
            return c
    # heuristic: longest median string column
    best, best_med = None, -1
    for c in df.columns:
        if df[c].dtype == object:
            med = df[c].astype(str).str.len().median()
            if med > best_med:
                best, best_med = c, med
    if best is None:
        raise ValueError("Could not detect sequence column in design TSV.")
    return best


def sort_candidates(df: pd.DataFrame) -> pd.DataFrame:
    rank_col = None
    for c in df.columns:
        if c.lower() in {"rank", "ranking"}:
            rank_col = c
            break
    if rank_col:
        return df.sort_values(rank_col, ascending=True).reset_index(drop=True)

    score_col = None
    for c in df.columns:
        if "score" in c.lower():
            score_col = c
            break
    if score_col:
        # assume higher is better
        return df.sort_values(score_col, ascending=False).reset_index(drop=True)

    return df.reset_index(drop=True)


def map_amplicon_seeded(
    amp_seq: str,
    target: str,
    k_list=(31, 25, 21, 17),
    diag_bin: int = 5,
    min_hits: int = 3,
    max_hits_per_kmer: int = 10,
):
    """
    Seed mapping onto TARGET gene sequence:
      - collect (qpos, tpos) for k-mer matches (forward and reverse complement)
      - choose dominant diagonal (tpos - qpos)
      - infer mapped span on target
    """
    amp = norm(amp_seq)
    if len(amp) < min(k_list):
        return None

    for orient, q in (("fwd", amp), ("rc", rc(amp))):
        for k in k_list:
            if k > len(q):
                continue

            hits = []  # (qpos, tpos)
            for qpos in range(0, len(q) - k + 1):
                kmer = q[qpos:qpos + k]
                if "N" in kmer:
                    continue
                if is_low_complexity(kmer):
                    continue

                start = 0
                found = 0
                while True:
                    tpos = target.find(kmer, start)
                    if tpos == -1:
                        break
                    hits.append((qpos, tpos))
                    found += 1
                    if found >= max_hits_per_kmer:
                        break
                    start = tpos + 1

            if len(hits) < min_hits:
                continue

            bins = Counter()
            for qpos, tpos in hits:
                diag = tpos - qpos
                b = int(round(diag / diag_bin))
                bins[b] += 1

            best_bin, best_ct = bins.most_common(1)[0]
            if best_ct < min_hits:
                continue

            diag_center = best_bin * diag_bin
            clustered = [(qpos, tpos) for (qpos, tpos) in hits if abs((tpos - qpos) - diag_center) <= diag_bin]
            if len(clustered) < min_hits:
                continue

            t_starts = [tpos for _, tpos in clustered]
            start = min(t_starts)
            end = max(t_starts) + k
            if start < 0 or end > len(target) or end <= start:
                continue

            return {
                "start": int(start),
                "end": int(end),
                "orient": orient,
                "k_used": int(k),
                "hits": int(len(clustered)),
                "diag_center": int(diag_center),
            }

    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--amplicons", required=True, help="kulkarni_amplicons.tsv")
    ap.add_argument("--results", required=True, help="results folder root (e.g., ./results)")
    ap.add_argument("--out_prefix", default="kulkarni_eval_v4_1", help="prefix for output TSVs")
    ap.add_argument("--skip_gene", action="append", default=[], help="Gene symbols to skip (repeatable)")
    ap.add_argument("--topN", type=int, default=10, help="Top-N candidates to annotate per run")
    ap.add_argument("--min_overlap_bp", type=int, default=50, help="Min overlap bp to count hit")
    ap.add_argument("--seed_k", type=int, default=31, help="Largest seed k (also tries 25,21,17)")
    args = ap.parse_args()

    amps = pd.read_csv(args.amplicons, sep="\t", dtype=str)
    need = {"Gene", "amplicon", "amp_seq", "publication_id", "Flybase_id"}
    missing = need - set(amps.columns)
    if missing:
        raise SystemExit(f"ERROR: missing columns in amplicons TSV: {sorted(missing)}")

    amps["amp_seq_norm"] = amps["amp_seq"].apply(norm)
    amps["amp_len"] = amps["amp_seq_norm"].str.len().astype(int)

    res_root = Path(args.results)
    if not res_root.exists():
        raise SystemExit(f"ERROR: results folder not found: {res_root}")

    skip = set(args.skip_gene)
    k_list = (args.seed_k, 25, 21, 17)

    per_run_rows = []
    top_rows = []

    for (gene, fbgn), sub in amps.groupby(["Gene", "Flybase_id"], dropna=False):
        if gene in skip:
            continue

        gene_dir = res_root / f"{gene}_{fbgn}"
        if not gene_dir.exists():
            continue

        # Amp labels from publication_id (Amp1/Amp2/Amp3)
        amp_map = {}
        for _, r in sub.iterrows():
            amp_label = str(r["publication_id"]).strip()
            amp_map[amp_label] = {
                "amplicon_id": r["amplicon"],
                "seq": r["amp_seq_norm"],
                "len": int(r["amp_len"]),
            }

        for run_dir in sorted([p for p in gene_dir.iterdir() if p.is_dir()]):
            tsv = find_design_tsv(run_dir)
            if not tsv:
                continue

            # IMPORTANT: load true target gene (other_files/target.fa)
            target = load_true_target_sequence(run_dir)
            if not target:
                continue

            # map each amplicon onto the target gene
            mapped = {amp_label: map_amplicon_seeded(info["seq"], target, k_list=k_list)
                      for amp_label, info in amp_map.items()}

            df = pd.read_csv(tsv, sep="\t", dtype=str)
            if df.empty:
                continue

            seq_col = detect_seq_col(df)
            df["cand_seq"] = df[seq_col].apply(norm)
            df = sort_candidates(df)
            df["rank"] = range(1, len(df) + 1)

            # map candidates to target by substring (fwd or rc)
            cand_starts, cand_ends = [], []
            for s in df["cand_seq"]:
                i = target.find(s)
                if i != -1:
                    cand_starts.append(i)
                    cand_ends.append(i + len(s))
                    continue
                s_rc = rc(s)
                j = target.find(s_rc)
                if j != -1:
                    cand_starts.append(j)
                    cand_ends.append(j + len(s_rc))
                    continue
                cand_starts.append(None)
                cand_ends.append(None)

            df["cand_start"] = cand_starts
            df["cand_end"] = cand_ends

            summary = {"Gene": gene, "Flybase_id": fbgn, "run_folder": run_dir.name, "design_tsv": str(tsv)}
            m = re.search(r"rnai_len(\d+)", run_dir.name)
            summary["rnai_length_run"] = int(m.group(1)) if m else None

            for amp_label in sorted(amp_map.keys()):
                mp = mapped.get(amp_label)
                if mp is None:
                    summary[f"{amp_label}_mapped"] = False
                    summary[f"{amp_label}_best_rank_overlap"] = None
                    summary[f"{amp_label}_in_topN_overlap"] = False
                    continue

                a_start, a_end = mp["start"], mp["end"]
                summary[f"{amp_label}_mapped"] = True
                summary[f"{amp_label}_map_start"] = a_start
                summary[f"{amp_label}_map_end"] = a_end
                summary[f"{amp_label}_map_orient"] = mp["orient"]
                summary[f"{amp_label}_map_k"] = mp["k_used"]
                summary[f"{amp_label}_map_hits"] = mp["hits"]

                hits = df[df["cand_start"].notna()].copy()
                if hits.empty:
                    summary[f"{amp_label}_best_rank_overlap"] = None
                    summary[f"{amp_label}_in_topN_overlap"] = False
                    continue

                hits["ov_bp"] = hits.apply(
                    lambda r: overlap(a_start, a_end, int(r["cand_start"]), int(r["cand_end"])),
                    axis=1
                )
                hits = hits[hits["ov_bp"] >= args.min_overlap_bp]
                best_rank = int(hits["rank"].min()) if not hits.empty else None
                summary[f"{amp_label}_best_rank_overlap"] = best_rank
                summary[f"{amp_label}_in_topN_overlap"] = bool((hits["rank"] <= args.topN).any())

            per_run_rows.append(summary)

            # annotate topN candidates with overlap bp to each amplicon
            top = df.head(args.topN).copy()
            for amp_label, mp in mapped.items():
                if mp is None:
                    top[f"ov_{amp_label}"] = 0
                    continue
                a_start, a_end = mp["start"], mp["end"]
                top[f"ov_{amp_label}"] = top.apply(
                    lambda r: 0 if pd.isna(r["cand_start"])
                    else overlap(a_start, a_end, int(r["cand_start"]), int(r["cand_end"])),
                    axis=1
                )
            top["Gene"] = gene
            top["Flybase_id"] = fbgn
            top["run_folder"] = run_dir.name
            top_rows.append(top)

    out1 = Path(f"{args.out_prefix}_per_run_summary.tsv")
    pd.DataFrame(per_run_rows).to_csv(out1, sep="\t", index=False)
    print(f"Wrote: {out1}")

    if top_rows:
        out2 = Path(f"{args.out_prefix}_top_candidates.tsv")
        pd.concat(top_rows, ignore_index=True).to_csv(out2, sep="\t", index=False)
        print(f"Wrote: {out2}")
    else:
        print("NOTE: No top candidate tables produced.")


if __name__ == "__main__":
    main()