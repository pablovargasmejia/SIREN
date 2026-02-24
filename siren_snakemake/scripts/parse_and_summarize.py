#!/usr/bin/env python3
import argparse, csv, os

def parse_rnahybrid_blocks(infile, sirna_size, min_align_length):
    if min_align_length is None:
        min_align_length = sirna_size - 4

    block = []
    inside = False
    for line in infile:
        if line.startswith("target:"):
            if block:
                aligned = block[11].replace(" ", "").strip()
                if len(aligned) >= min_align_length:
                    yield block
            block = []
            inside = True
        if inside:
            block.append(line)

    if block:
        aligned = block[11].replace(" ", "").strip()
        if len(aligned) >= min_align_length:
            yield block

def summarize_offtargets(filtered_txt, tsv_out):
    off = {}
    target = None
    with open(filtered_txt) as fh:
        for raw in fh:
            line = raw.strip()
            if line.startswith("target:"):
                target = line.split(": ", 1)[1].strip()
                off.setdefault(target, {"count":0, "sirnas":[]})
            elif line.startswith("miRNA :") and target:
                sirna = line.split(": ", 1)[1].strip()
                off[target]["count"] += 1
                off[target]["sirnas"].append(sirna)
            elif line == "":
                target = None

    items = sorted(off.items(), key=lambda x: x[1]["count"], reverse=True)
    with open(tsv_out, "w", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(["Off target","siRNA number","siRNA names"])
        for t, d in items:
            w.writerow([t, d["count"], ", ".join(d["sirnas"])])

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--all_targets", required=True)
    ap.add_argument("--sirna_size", type=int, default=21)
    ap.add_argument("--min_align_length", type=int, default=None)
    ap.add_argument("--off_targets_results", required=True)
    ap.add_argument("--off_targets_summary", required=True)
    args = ap.parse_args()

    os.makedirs(os.path.dirname(args.off_targets_results), exist_ok=True)

    with open(args.all_targets) as inp, open(args.off_targets_results, "w") as out:
        for block in parse_rnahybrid_blocks(inp, args.sirna_size, args.min_align_length):
            out.writelines(block)

    summarize_offtargets(args.off_targets_results, args.off_targets_summary)

if __name__ == "__main__":
    main()
