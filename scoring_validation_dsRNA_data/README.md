# SIREN scoring system validation against a dsRNA off‑target dataset (Drosophila; Kulkarni 2006 study)

## Context and question

To evaluate whether **SIREN** tends to avoid designing dsRNAs that resemble experimentally observed off‑target–prone reagents, I benchmarked SIREN’s ranked dsRNA candidates against a Drosophila melanogaster cell‑based dataset reported in *Evidence of off‑target effects associated with long dsRNAs in Drosophila melanogaster cell‑based assays* (Kulkarni et al.). For each tested gene, that study provides three long-dsRNA amplicons (labeled **Amp1**, **Amp2**, and **Amp3**). In the way this dataset is used here, **Amp1** is treated as the “problematic” amplicon (associated with off‑target effects in the study), whereas **Amp2/Amp3** are treated as the comparatively “clean” regions for that gene.

The central question is simple: when SIREN proposes its best dsRNA designs for the same gene, does it avoid the Amp1 region and (when it overlaps the published amplicon regions at all) does it prefer Amp2/Amp3?

## What was evaluated

This benchmark covers six genes with completed runs (**CG12155**, **CG30421**, **CG32791**, **CG3563**, **Smox**, and **trio**). The gene **shot** was excluded at the time of evaluation because its SIREN run was still in progress and was expected to be unusually large.

Each gene was run three times in SIREN, one run per experimental amplicon length. In the results folder, each run corresponds to a distinct output directory whose name includes the relevant `rnai_len<LEN>` tag; the analysis below therefore evaluates results **per run**, and automatically keeps results separated by dsRNA length.

## How overlap between SIREN designs and Amp1/Amp2/Amp3 was assessed

For each finished run, the “true” target gene sequence used by SIREN was taken from `<run_outdir>/other_files/target.fa` (not from the off‑target database, which is stored elsewhere). Each published amplicon sequence (Amp1/Amp2/Amp3) was then mapped onto that target sequence using a seed‑based mapping strategy (progressively trying k‑mers 31→25→21→17 while filtering low‑complexity seeds). This step yields target‑coordinate intervals for Amp1/Amp2/Amp3 on the target sequence that SIREN actually used.

Next, each ranked SIREN candidate dsRNA sequence was mapped back onto the same target sequence by exact substring match (forward strand or reverse‑complement). A candidate was considered to “hit” an amplicon region if its mapped interval overlapped the mapped Amp interval by **at least 50 bp** (the analysis threshold used here). With this definition, the output statistics track (i) whether the **#1** ranked design overlaps Amp1/Amp2/Amp3 and (ii) whether any of the **top 20** designs overlap Amp1/Amp2/Amp3, as well as the best (lowest) rank at which an overlap occurs.

## Overall results

Across the **18 finished runs** (six genes × three lengths), SIREN’s **top‑ranked** design did **not** overlap Amp1 in any run under the ≥50 bp criterion. When the top design overlapped any published amplicon region, it overlapped Amp2 or Amp3 instead.

|   Finished runs |   Top-1 overlaps Amp1 |   Top-1 overlaps Amp2 |   Top-1 overlaps Amp3 |   Top-1 overlaps any amplicon |   Top-20 contains Amp1 |   Top-20 contains Amp2 |   Top-20 contains Amp3 |
|----------------:|----------------------:|----------------------:|----------------------:|------------------------------:|-----------------------:|-----------------------:|-----------------------:|
|              18 |                     0 |                     2 |                     5 |                             7 |                      1 |                      3 |                      5 |

## Gene‑level interpretation

For **CG12155 (FBgn0029957)**, SIREN’s #1 design overlapped **Amp3** in the two longer runs (rnai_len 516 and 491), while Amp1 only appeared far down the ranking (best overlap ranks 33 and 34, respectively). The short run (rnai_len 202) did not overlap any of the three published amplicons at ≥50 bp within the top 20 candidates, suggesting SIREN’s top designs in that run sit outside the tested amplicon regions under the current overlap threshold.

For **CG32791 (FBgn0052791)**, SIREN’s #1 design overlapped **Amp2** for rnai_len 455 and 353. In the longest run (rnai_len 509), the #1 design did not overlap any published amplicon region, but Amp2 still appeared near the top of the ranking (best overlap rank 3). This pattern is consistent with SIREN preferentially selecting the “clean” region Amp2 when it selects among the published regions, and otherwise selecting a different region of the transcript.

For **Smox (FBgn0025800)**, SIREN’s #1 design overlapped **Amp3** in all three runs (rnai_len 506, 342, and 362). Notably, the only instance where Amp1 appeared in the top 20 candidates across the entire benchmark occurred in the Smox rnai_len 506 run, where Amp1 overlapped at rank 7 while Amp3 overlapped at rank 1. Even in this worst case, Amp1 was not the top choice.

For **CG30421 (FBgn0050421)**, **CG3563 (FBgn0263929)**, and **trio (FBgn0024277)**, none of the top 20 candidates overlapped Amp1/Amp2/Amp3 by ≥50 bp in any run. In practical terms, these genes are not informative for the specific “Amp1 versus Amp2/Amp3 preference” question under the current overlap rule, because SIREN’s best designs consistently fall outside the experimentally tested amplicon regions.

The table below summarizes these patterns compactly.

| Gene    | FlyBase ID   |   Runs |   Top1 Amp1 |   Top1 Amp2 |   Top1 Amp3 |   Top20 Amp1 |   Top20 Amp2 |   Top20 Amp3 | Best rank Amp1   | Best rank Amp2   | Best rank Amp3   |
|:--------|:-------------|-------:|------------:|------------:|------------:|-------------:|-------------:|-------------:|:-----------------|:-----------------|:-----------------|
| CG12155 | FBgn0029957  |      3 |           0 |           0 |           2 |            0 |            0 |            2 | 33.0             |                  | 1.0              |
| CG30421 | FBgn0050421  |      3 |           0 |           0 |           0 |            0 |            0 |            0 |                  |                  |                  |
| CG32791 | FBgn0052791  |      3 |           0 |           2 |           0 |            0 |            3 |            0 |                  | 1.0              |                  |
| CG3563  | FBgn0263929  |      3 |           0 |           0 |           0 |            0 |            0 |            0 |                  |                  |                  |
| Smox    | FBgn0025800  |      3 |           0 |           0 |           3 |            1 |            0 |            3 | 7.0              |                  | 1.0              |
| trio    | FBgn0024277  |      3 |           0 |           0 |           0 |            0 |            0 |            0 |                  |                  |                  |

## Conclusions

Under a conservative overlap definition (≥50 bp overlap on the SIREN target sequence), SIREN’s ranked dsRNA designs show a strong qualitative agreement with the dataset’s intended benchmark logic: **Amp1 is not selected as the top design in any finished run**, and when SIREN’s top design overlaps a published amplicon region, it overlaps **Amp2 or Amp3** rather than Amp1. The primary limitation is that, for half of the evaluated genes, SIREN selects top candidates outside all three published amplicon regions, making those genes **inconclusive** for distinguishing Amp1 from Amp2/Amp3 using overlap alone (not necessarily indicative of poor design choices).

## Notes and limitations

These results quantify positional overlap with experimentally tested amplicon regions; they do not directly model phenotype. The strict ≥50 bp overlap threshold is intentionally conservative and could be reduced to capture partial overlaps, at the cost of potentially counting weak or incidental overlaps. Finally, the **shot** gene is excluded here until its run finishes; incorporating it later may be informative, but it is expected to be computationally heavy due to transcript complexity and size.

## References
Kulkarni, M.M., Booker, M., Silver, S.J., Friedman, A., Hong, P., Perrimon, N., et al. (2006) Evidence of off-target effects associated with long dsRNAs in Drosophila melanogaster cell-based assays. Nature Methods, 3, 833–838. https://doi.org/10.1038/nmeth935.
