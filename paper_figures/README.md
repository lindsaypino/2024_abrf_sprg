# paper_figures/

Canonical, self-contained code for each paper figure. Each analysis figure is one
notebook that reads the raw participant returns directly and writes its own outputs —
no shared parsing step, no `H:/My Drive` paths. Not every figure is code-driven (some
are schematics or instrument screenshots); those are listed here with their source.

## Figure index

| Figure | Title | Source | Status |
|--------|-------|--------|--------|
| Fig 1 | Sample / experiment design | schematic (drawn) — no notebook | placeholder |
| **Fig 2** | Single-dataset example of the Fig 3 analysis (site 42 / Talus) | [`fig2_single_dataset_example.ipynb`](fig2_single_dataset_example.ipynb) → **FIG02** | ✅ done (DIA-NN `report.pr_matrix`) |
| **Fig 3** | Expected vs. observed (recalculated from peptide quant) | [`fig3_expected_vs_observed.ipynb`](fig3_expected_vs_observed.ipynb) → **FIG03A** facets + **FIG03B** recalculated box-whisker | ✅ done |
| **Fig 4** | Observed log2(A/C) density by species | [`fig4_ratio_density_by_species.ipynb`](fig4_ratio_density_by_species.ipynb) → **FIG04** | ✅ done |
| **Fig 5** | Human/bovine (vs. trout) homology challenge | [`fig5_homology.ipynb`](fig5_homology.ipynb) → **FIG05** (KDE density of log2(A/C) by homology class, DIA-NN). Companion Skyline SAAV chromatograms = manual instrument screenshot. | ✅ code done (DIA-NN) |
| Fig 6? | Reanalysis | TBD | placeholder |
| **Supp** | Summary of participant-reported sample ratios | [`supp_reported_sample_ratios.ipynb`](supp_reported_sample_ratios.ipynb) → **SUPP_reported_ratios** | ✅ done |

Figure numbering is still in flux between the outline and the placeholder deck; names
here are descriptive.

## Notebooks

- **`fig2_single_dataset_example.ipynb`** — single-dataset LFQbench-style worked example
  (observed log2 fold-change vs. abundance + per-species boxplots, faceted by A/B, A/C, B/C)
  for the site-42 run → `FIG02_single_dataset_example`. Species are assigned by in-silico
  tryptic digest of `Combined_proteomes.fasta` (unique peptides only; map cached to `data/`).
  Reads the DIA-NN `report.pr_matrix.tsv` (precursor-level, collapsed to peptide by summing
  across charges). `INPUT`/`COLS` at the top parameterize the input if it changes again.
- **`fig3_expected_vs_observed.ipynb`** — from each participant's `*-pepQuant` peptide
  intensities. Produces two figures:
  - `FIG03A_expected_vs_observed_facets` — per-submission observed log2(A/C) vs. abundance.
  - `FIG03B_recalculated_ratios` — A/B, B/C, A/C box-and-whisker by species (median of
    per-peptide ratios, datatype-aware; robust to high-abundance / human-cow-shared
    peptides that bias a ratio-of-summed-intensities).
- **`fig4_ratio_density_by_species.ipynb`** — per-species KDE of observed log2(A/C),
  one curve per submission → `FIG04_ratio_density_by_species`.
- **`fig5_homology.ipynb`** — the human/bovine homology challenge → `FIG05_homology_challenge`.
  KDE density of observed log2(A/C) coloured by peptide homology class. Human-unique peaks at
  +3.9 and bovine-unique at −3.2 (expected), but human+bovine-**shared** peptides collapse to
  a broad peak at ≈0 (ratio 1) — because human+cow sum to a constant 50 % per sample they
  can't resolve the two species, biasing protein roll-up. Homology class is sequence-level
  (in-silico tryptic digest of the combined FASTA, cached to `data/`); DIA-NN's `Protein.Names`
  is deliberately **not** used for this (it collapses shared peptides onto one protein group
  and mislabels ~2/3 of them as unique). Skyline SAAV chromatograms added manually.
- **`fig5_homology_saav_helper.ipynb`** — helper (no figure): lists candidate single-amino-
  acid-variant peptide pairs for the Skyline chromatogram panel. Finds human-unique/
  bovine-unique peptides differing by exactly one residue (position-wildcard index) that
  co-elute within ~30 s (median RT from `report.tsv`), flags pairs where the human peptide's
  A/C is flattened toward 1 (co-eluting-SAAV interference), and writes
  `output/fig5_saav_candidates.csv`. Reimplements the neighbor search in
  `bin/06_homology/peptide_ratios.ipynb`, adding the RT-proximity filter.
- **`supp_reported_sample_ratios.ipynb`** — the ratios each participant *reported* in
  their `*-sampleRatios` file → `SUPP_reported_ratios` (companion to FIG03B).

All carry the per-submission **curation** inline (verified data types; fake-site
exclusions 05/42/Templates; correct file choices for 02/03/14; A/B splits for 07/10/14;
site-10-QEx and site-09 column fixups). See each notebook's header cell.

## Layout

```
paper_figures/
  fig2_single_dataset_example.ipynb      -> output/FIG02_*
  fig3_expected_vs_observed.ipynb        -> output/FIG03A_*, output/FIG03B_*
  fig4_ratio_density_by_species.ipynb    -> output/FIG04_*
  fig5_homology.ipynb                    -> output/FIG05_*
  fig5_homology_saav_helper.ipynb        -> output/fig5_saav_candidates.csv (Skyline picks)
  supp_reported_sample_ratios.ipynb      -> output/SUPP_*
  output/     # rendered PNG + PDF (committed)
  data/       # tidy long-format CSVs (regenerated by the notebooks; git-ignored)
```

## Running

Data source (machine-specific): `D:/2022 Multi-Species Standard Study/`.
Open a notebook in Jupyter and Run-All, or headless:

```
python -m nbconvert --to notebook --execute --inplace paper_figures/fig3_expected_vs_observed.ipynb
```

Requires: pandas, numpy, matplotlib, scipy, openpyxl.
