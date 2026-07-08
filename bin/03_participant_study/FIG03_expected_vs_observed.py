"""
FIG03: Expected vs. observed ratios, per participant submission.

Observed log2(A/C) per peptide plotted against abundance, with a dashed line at the
design-expected log2(A/C) for each species. Two views:

  FIG03A  per-submission facets : x = abundance (type-aware), y = log2(A/C) observed,
                                  points colored by species, dashed = expected per species.
  FIG03B  per-species KDE       : density of observed log2(A/C), one curve per submission,
                                  vertical dashed = expected.

Built to mirror bin/04_summary/FIG04_reported_vs_recalculated.py — SAME curated file map,
SAME per-file fixups, SAME fake-data exclusions. Each submission (incl. multi-workflow
A/B variants) is its own entry.

DATA-TYPE HANDLING (verified empirically from the quantity magnitudes, not assumed):
  raw    {01,02,03,04,07A,07B,09,11,12,14A,14B} : intensities  -> x=log2(B), y=log2(A/C)
  log2   {06}   : values already in log2 space   -> x=B,           y=A-C
  scaled {10A,10B} : normalised to 300 % (per-peptide A+B+C=300) -> y=log2(A/C) is valid,
                     but x=log2(B) is a percentage, NOT an intensity -> abundance x flagged
                     non-comparable (kept in ratio/KDE panels, marked in facet title).
  counts {41}   : spectral-count-like integers    -> x=log2(B+1), y=log2((A+1)/(C+1));
                     ratio is count-based -> marked non-comparable in title/caption.

EXCLUDED: sites 05 and 42 (only 'FAKEPEPTIDE' placeholder pepQuant) and 14/Templates
(fake copy) — never enter the curated map.
"""

import os
import re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.stats import gaussian_kde

DATA_DIR = r"D:/2022 Multi-Species Standard Study"
FIGURES  = r"D:/2024_abrf_sprg/figures"
RESULTS  = r"D:/2024_abrf_sprg/results"
os.makedirs(FIGURES, exist_ok=True)
os.makedirs(RESULTS, exist_ok=True)

# label -> pepQuant file (relative to DATA_DIR). Identical curation to FIG04.
CURATED = {
    "01":  "01/01-pepQuant.tsv",
    "02":  "02/02-pepQuant.csv",
    "03":  "03/03-pepQuant.tsv",
    "04":  "04/04-pepQuant.tsv",
    "06":  "06/06-pepQuant.tsv",
    "07A": "07/07-pepQuant.tsv",
    "07B": "07/07-pepQuant_ShortGrad.tsv",
    "09":  "09/09-pepQuant.tsv",
    "10A": "10/10-QEx_pepQuant.tsv",
    "10B": "10/10-Fusion_pepQuant.tsv",
    "11":  "11/11-pepQuant.tsv",
    "12":  "12/12-pepQuant.tsv",
    "14A": "14/14_pepQuant_sPRG_lumos_DIA_1x8mzStag_3x4mzGPFlibrary.txt",
    "14B": "14/14_pepQuant_sPRG_lumos_DIA_2x4mzStag_Prosit_library.txt",
    "41":  "41/41-pepQuant.tsv",
}
ORDER = list(CURATED.keys())

# verified per-submission data type (see module docstring)
DATATYPE = {k: "raw" for k in CURATED}
DATATYPE["06"] = "log2"
DATATYPE["10A"] = "scaled"
DATATYPE["10B"] = "scaled"
DATATYPE["41"] = "counts"
NONCOMPARABLE_X = {"scaled", "counts"}   # abundance axis not comparable across sites

SPECIES_MAP = {
    "cow": "Bovine", "bovin": "Bovine", "bovine": "Bovine", "bos taurus": "Bovine",
    "human": "Human", "homo sapiens": "Human",
    "trout": "Trout", "salvelinus namaycush": "Trout", "salnm": "Trout",
}
SPECIES_ORDER = ["Bovine", "Human", "Trout"]
COLORS = {"Bovine": "#E69F00", "Human": "#56B4E9", "Trout": "#009E73"}

# design mix (peptide-level %): expected log2(A/C)
MIX = {"A": {"Trout": 50, "Human": 45, "Bovine": 5},
       "C": {"Trout": 50, "Human": 3,  "Bovine": 47}}
EXPECTED_AC = {sp: np.log2(MIX["A"][sp] / MIX["C"][sp]) for sp in SPECIES_ORDER}


def _read(path):
    if path.lower().endswith((".xlsx", ".xls")):
        return pd.read_excel(path)
    sep = "," if path.lower().endswith(".csv") else "\t"
    return pd.read_csv(path, sep=sep, engine="python", on_bad_lines="skip")


def load_peptides():
    """One tidy frame: label, datatype, species, x_metric, y_metric (log2 A/C)."""
    frames = []
    for label, pq in CURATED.items():
        df = _read(os.path.join(DATA_DIR, pq))
        df.columns = [str(c).strip() for c in df.columns]

        # FIXUP: site 10-QEx 3rd quantity col mislabeled "Sample B Quantity" (dup) = Sample C
        if "Sample B Quantity.1" in df.columns and "Sample C Quantity" not in df.columns:
            df = df.rename(columns={"Sample B Quantity.1": "Sample C Quantity"})

        scol = "Species" if "Species" in df.columns else \
               ("PG.Organisms" if "PG.Organisms" in df.columns else None)
        if scol is None:
            raise ValueError(f"{label}: no species column ({list(df.columns)})")
        species = df[scol].map(lambda x: SPECIES_MAP.get(str(x).strip().lower()))

        # per-sample quantity, averaging any replicate columns (site 09 A_R1/R2/R3 ...)
        q = {}
        for L in ("A", "B", "C"):
            cols = [c for c in df.columns if re.search(rf"Sample {L}(_R\d+)? Quantity", c)]
            if not cols:
                raise ValueError(f"{label}: no 'Sample {L} Quantity' column in {pq}")
            q[L] = df[cols].apply(pd.to_numeric, errors="coerce").mean(axis=1)
        qA, qB, qC = q["A"], q["B"], q["C"]

        dt = DATATYPE[label]
        if dt == "raw" or dt == "scaled":
            keep = (qA > 0) & (qB > 0) & (qC > 0)
            x = np.log2(qB.where(keep))
            y = np.log2((qA / qC).where(keep))
        elif dt == "log2":                       # values already log2
            keep = qA.notna() & qB.notna() & qC.notna()
            x = qB.where(keep)
            y = (qA - qC).where(keep)
        elif dt == "counts":
            keep = qA.notna() & qB.notna() & qC.notna()
            x = np.log2(qB.where(keep) + 1)
            y = np.log2((qA.where(keep) + 1) / (qC.where(keep) + 1))

        sub = pd.DataFrame({"label": label, "datatype": dt, "species": species,
                            "x_metric": x, "y_metric": y})
        sub = sub[sub["species"].isin(SPECIES_ORDER)].dropna(subset=["x_metric", "y_metric"])
        frames.append(sub)
    return pd.concat(frames, ignore_index=True)


def fig_facets(p, out_png):
    n = len(ORDER)
    ncol, nrow = 4, int(np.ceil(n / 4))
    y_lo, y_hi = np.nanpercentile(p["y_metric"], [1, 99])
    ypad = 0.08 * (y_hi - y_lo)
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.6 * ncol, 2.8 * nrow),
                             squeeze=False, sharey=True)
    for i, label in enumerate(ORDER):
        ax = axes[divmod(i, ncol)]
        sub = p[p["label"] == label]
        dt = DATATYPE[label]
        for sp in SPECIES_ORDER:
            ss = sub[sub["species"] == sp]
            if not ss.empty:
                ax.scatter(ss["x_metric"], ss["y_metric"], s=6, alpha=0.35,
                           color=COLORS[sp], linewidth=0)
            ax.axhline(EXPECTED_AC[sp], ls="--", lw=1.0, color=COLORS[sp], zorder=3)
        tag = "" if dt == "raw" else f"  [{dt}]"
        ax.set_title(f"{label}{tag}", fontsize=10,
                     color=("black" if dt == "raw" else "#B00000"))
        ax.set_ylim(y_lo - ypad, y_hi + ypad)
        ax.grid(True, lw=0.3)
        if dt in NONCOMPARABLE_X:                       # abundance not comparable
            ax.set_facecolor("#f7f2ea")
    for j in range(n, nrow * ncol):
        axes[divmod(j, ncol)].axis("off")
    handles = [Line2D([0], [0], marker="o", ls="", color=COLORS[s], label=s) for s in SPECIES_ORDER]
    handles += [Line2D([0], [0], ls="--", color="0.3", label="expected log2(A/C)")]
    axes[divmod(n, ncol)].legend(handles=handles, loc="center", frameon=False, fontsize=10) \
        if n < nrow * ncol else fig.legend(handles=handles, loc="lower center", ncol=4)
    fig.supxlabel("abundance  (raw/scaled: log2 B  |  log2: B  |  counts: log2(B+1))   "
                  "— shaded panels = abundance not comparable")
    fig.supylabel("observed log2(A/C)   (log2 site: A-C;  counts: log2((A+1)/(C+1)))")
    fig.suptitle("FIG03  Expected vs. observed log2(A/C) per submission", fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.98))
    fig.savefig(out_png, dpi=200); fig.savefig(out_png.replace(".png", ".pdf"))
    plt.close(fig)
    print(f"  wrote {out_png}")


def fig_kde(p, out_png, min_pts=50):
    vals = p["y_metric"].replace([np.inf, -np.inf], np.nan).dropna()
    x_lo, x_hi = np.nanpercentile(vals, [0.5, 99.5])
    pad = 0.05 * (x_hi - x_lo)
    grid = np.linspace(x_lo - pad, x_hi + pad, 512)
    cmap = plt.get_cmap("tab20")
    colmap = {lab: cmap(i % 20) for i, lab in enumerate(ORDER)}

    fig = plt.figure(figsize=(13, 8))
    gs = fig.add_gridspec(2, 2)
    axes = {"Bovine": fig.add_subplot(gs[0, 0]), "Human": fig.add_subplot(gs[0, 1]),
            "Trout": fig.add_subplot(gs[1, 0])}
    ax_leg = fig.add_subplot(gs[1, 1]); ax_leg.axis("off")

    for sp, ax in axes.items():
        ax.axvline(EXPECTED_AC[sp], ls="--", lw=1.2, color="black", zorder=5)
        for lab in ORDER:
            v = p[(p["species"] == sp) & (p["label"] == lab)]["y_metric"]
            v = v.replace([np.inf, -np.inf], np.nan).dropna().values
            if len(v) < min_pts:
                continue
            dens = gaussian_kde(v)(grid)
            ax.plot(grid, dens, lw=2.4, alpha=0.8, color=colmap[lab],
                    ls=("--" if DATATYPE[lab] == "counts" else "-"))
        ax.set_title(sp, color=COLORS[sp]); ax.grid(True, lw=0.3)
        ax.set_xlim(grid[0], grid[-1]); ax.set_ylabel("density")
        ax.set_xlabel("observed log2(A/C)")
    handles = [Line2D([0], [0], color=colmap[l], lw=3,
                      label=l + (" *" if DATATYPE[l] in ("counts",) else "")) for l in ORDER]
    handles += [Line2D([0], [0], ls="--", color="black", lw=1.2, label="expected log2(A/C)"),
                Line2D([0], [0], ls="--", color="0.4", lw=2.4, label="* count data (dashed)")]
    ax_leg.legend(handles=handles, loc="center", ncol=2, frameon=False, title="submission")
    fig.suptitle("FIG03  Observed log2(A/C) density by species — dashed = design-expected", fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(out_png, dpi=200); fig.savefig(out_png.replace(".png", ".pdf"))
    plt.close(fig)
    print(f"  wrote {out_png}")


if __name__ == "__main__":
    p = load_peptides()
    p.to_csv(os.path.join(RESULTS, "fig03_peptide_ratios_long.csv"), index=False)
    print("expected log2(A/C):", {k: round(v, 3) for k, v in EXPECTED_AC.items()})
    print(f"{len(p)} peptide points, {p['label'].nunique()} submissions -> {sorted(p['label'].unique())}")
    print(p.groupby("label").size().to_string())
    fig_facets(p, os.path.join(FIGURES, "FIG03A_expected_vs_observed_facets.png"))
    fig_kde(p, os.path.join(FIGURES, "FIG03B_ratio_density_by_species.png"))
