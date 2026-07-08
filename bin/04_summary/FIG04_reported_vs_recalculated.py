"""
FIG04 (SUPP): Summary of participant ratios — reported vs. recalculated.

Two box-and-whisker figures over the 3 pairwise sample comparisons (A/B, B/C, A/C),
split by species (Bovine, Human, Trout), with dashed lines at the design-expected ratios:

  FIG_A  reported_ratios      : the ratios each participant reported (their *-sampleRatios files)
  FIG_B  recalculated_ratios  : ratios we recompute from each participant's peptide quantities
                                (sum of peptide intensity per species per sample -> A/B, B/C, A/C)

CURATION NOTES (examined every file in each subdirectory; see conversation log):
  * Sites 05 and 42 are EXCLUDED entirely: both submitted only placeholder template data
    ("FAKEPEPTIDE" in pepQuant; "*this is fake data" in sampleRatios). No usable quant.
  * Site 02 : the .tsv sampleRatios is the fake template -> we use 02-sampleRatios.xlsx (real).
  * Site 03 : the plain .tsv is the fake template -> we use 03-sampleRatios_completed.tsv.
  * Site 14 : Templates/14-sampleRatios.tsv is the fake template -> excluded; the two real
              workflow variants are kept as 14A / 14B.
  * Multi-workflow labs are kept as separate entries (A/B suffix), matching Table 1:
        07A = 120 min (long) gradient      07B = 90 min (ShortGrad)
        10A = Q Exactive ("System 1")      10B = Orbitrap Fusion + FAIMS ("System 2")
        14A = 8 m/z staggered + GPF library  14B = 4 m/z staggered + Prosit library
"""

import os
import re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

DATA_DIR = r"D:/2022 Multi-Species Standard Study"
RESULTS  = r"D:/2024_abrf_sprg/results"
FIGURES  = r"D:/2024_abrf_sprg/figures"
os.makedirs(RESULTS, exist_ok=True)
os.makedirs(FIGURES, exist_ok=True)

# label -> curated (sampleRatios file, pepQuant file). Paths relative to DATA_DIR.
CURATED = {
    "01":  ("01/01-sampleRatios.tsv",                                          "01/01-pepQuant.tsv"),
    "02":  ("02/02-sampleRatios.xlsx",                                         "02/02-pepQuant.csv"),   # .tsv is fake
    "03":  ("03/03-sampleRatios_completed.tsv",                                "03/03-pepQuant.tsv"),   # plain .tsv is fake
    "04":  ("04/04-sampleRatios.tsv",                                          "04/04-pepQuant.tsv"),
    "06":  ("06/06-sampleRatios.tsv",                                          "06/06-pepQuant.tsv"),
    "07A": ("07/07-sampleRatios.tsv",                                          "07/07-pepQuant.tsv"),
    "07B": ("07/07-sampleRatios_ShortGrad.tsv",                                "07/07-pepQuant_ShortGrad.tsv"),
    "09":  ("09/09-sampleRatios.tsv",                                          "09/09-pepQuant.tsv"),
    "10A": ("10/10-QEx_sampleRatios.tsv",                                      "10/10-QEx_pepQuant.tsv"),
    "10B": ("10/10-Fusion_sampleRatios.tsv",                                   "10/10-Fusion_pepQuant.tsv"),
    "11":  ("11/11-sampleRatios.tsv",                                          "11/11-pepQuant.tsv"),
    "12":  ("12/12-sampleRatios.tsv",                                          "12/12-pepQuant.tsv"),
    "14A": ("14/14-sampleRatios_sPRG_lumos_DIA_1x8mzStag_3x4mzGPFlibrary.txt", "14/14_pepQuant_sPRG_lumos_DIA_1x8mzStag_3x4mzGPFlibrary.txt"),
    "14B": ("14/14-sampleRatios_sPRG_lumos_DIA_2x4mzStag_Prosit_library.txt",  "14/14_pepQuant_sPRG_lumos_DIA_2x4mzStag_Prosit_library.txt"),
    "41":  ("41/41-sampleRatios.tsv",                                          "41/41-pepQuant.tsv"),
}

SPECIES_MAP = {  # normalize any spelling -> plotting label
    "cow": "Bovine", "bovin": "Bovine", "bovine": "Bovine", "bos taurus": "Bovine",
    "human": "Human", "homo sapiens": "Human",
    "trout": "Trout", "salvelinus namaycush": "Trout", "salnm": "Trout",
}
COMPARISONS = ["A/B", "B/C", "A/C"]
SPECIES_ORDER = ["Bovine", "Human", "Trout"]
COLORS = {"Bovine": "#E69F00", "Human": "#56B4E9", "Trout": "#009E73"}

EXPECTED = {  # design (peptide-level %): A=45/5/50, B=20/30/50, C=3/47/50 (human/cow/trout)
    "A/B": {"Bovine": 5/30,  "Human": 45/20, "Trout": 1.0},
    "B/C": {"Bovine": 30/47, "Human": 20/3,  "Trout": 1.0},
    "A/C": {"Bovine": 5/47,  "Human": 45/3,  "Trout": 1.0},
}


def _read(path):
    """Robust tabular read for pepQuant (skips the rare malformed row)."""
    if path.lower().endswith((".xlsx", ".xls")):
        return pd.read_excel(path)
    sep = "," if path.lower().endswith(".csv") else "\t"
    return pd.read_csv(path, sep=sep, engine="python", on_bad_lines="skip")


def _read_ratios(path):
    """sampleRatios reader: keep only the first 4 columns (Sample, A/B, B/C, A/C).
    Tolerates trailing tabs / trailing note columns / ragged rows."""
    if path.lower().endswith((".xlsx", ".xls")):
        df = pd.read_excel(path)
        df = df.iloc[:, :4]
        df.columns = ["Sample", "A/B", "B/C", "A/C"]
        return df
    recs = []
    with open(path, encoding="utf-8", errors="ignore") as fh:
        lines = [ln.rstrip("\n") for ln in fh if ln.strip()]
    for ln in lines[1:]:
        f = (ln.split("\t") + ["", "", "", ""])[:4]
        recs.append(f)
    return pd.DataFrame(recs, columns=["Sample", "A/B", "B/C", "A/C"])


def norm_species(s):
    return str(s).strip().lower()


# ---------------------------------------------------------------- reported
def load_reported():
    rows = []
    for label, (sr, _) in CURATED.items():
        df = _read_ratios(os.path.join(DATA_DIR, sr))
        sample_col = df.columns[0]  # first col holds the species name
        df["species"] = df[sample_col].map(lambda x: SPECIES_MAP.get(norm_species(x)))
        df = df.dropna(subset=["species"])
        for comp in COMPARISONS:
            if comp not in df.columns:
                continue
            for _, r in df.iterrows():
                val = pd.to_numeric(r[comp], errors="coerce")
                if pd.notna(val):
                    rows.append(dict(label=label, species=r["species"],
                                     comparison=comp, ratio=float(val)))
    return pd.DataFrame(rows)


# ------------------------------------------------------------ recalculated
def load_recalculated():
    rows = []
    for label, (_, pq) in CURATED.items():
        df = _read(os.path.join(DATA_DIR, pq))
        df.columns = [str(c).strip() for c in df.columns]

        # --- per-file column fixups (participant labeling errors) ---
        # Site 10-QEx: 3rd quantity column is mislabeled "Sample B Quantity" (dup);
        # pandas renames the duplicate to "Sample B Quantity.1" -> it is really Sample C.
        if "Sample B Quantity.1" in df.columns and "Sample C Quantity" not in df.columns:
            df = df.rename(columns={"Sample B Quantity.1": "Sample C Quantity"})
        # Site 09 mislabels "Peptide" as "Pepide"; site 14 uses "Sequence" (peptide id
        # is not needed for species-summed ratios, so no rename required here).

        # species column (site 07 uses PG.Organisms)
        scol = "Species" if "Species" in df.columns else \
               ("PG.Organisms" if "PG.Organisms" in df.columns else None)
        if scol is None:
            raise ValueError(f"{label}: no species column in {pq} ({list(df.columns)})")
        df["species"] = df[scol].map(norm_species).map(SPECIES_MAP)

        # per-sample quantity: average any replicate columns (site 09 has A_R1/R2/R3 ...)
        tot = {}
        for L in ("A", "B", "C"):
            qcols = [c for c in df.columns if re.search(rf"Sample {L}(_R\d+)? Quantity", c)]
            if not qcols:
                raise ValueError(f"{label}: no 'Sample {L} Quantity' column in {pq}")
            df[f"q{L}"] = df[qcols].apply(pd.to_numeric, errors="coerce").mean(axis=1)

        # sum peptide intensity per species, then take pairwise ratios
        g = df.dropna(subset=["species"]).groupby("species")[["qA", "qB", "qC"]].sum()
        for sp, r in g.iterrows():
            ratios = {"A/B": r.qA / r.qB, "B/C": r.qB / r.qC, "A/C": r.qA / r.qC}
            for comp, val in ratios.items():
                if np.isfinite(val) and val > 0:
                    rows.append(dict(label=label, species=sp, comparison=comp, ratio=val))
    return pd.DataFrame(rows)


# ------------------------------------------------------------------- plot
def boxwhisker(p, title, out_png):
    fig, axes = plt.subplots(3, 1, figsize=(5.4, 7.8), sharex=True)
    x_lim = (0.02, 30)
    for ax, comp in zip(axes, COMPARISONS):
        subc = p[p["comparison"] == comp]
        data = [subc.loc[subc["species"] == sp, "ratio"].values for sp in SPECIES_ORDER]
        ax.boxplot(data, orientation="horizontal", tick_labels=SPECIES_ORDER, widths=0.6, showfliers=False)
        ax.set_xscale("log")
        rng = np.random.default_rng(0)
        for i, sp in enumerate(SPECIES_ORDER, start=1):
            ss = subc[subc["species"] == sp]
            if not ss.empty:
                y = np.full(len(ss), i) + rng.uniform(-0.09, 0.09, size=len(ss))
                ax.scatter(ss["ratio"], y, s=18, alpha=0.8, color=COLORS[sp],
                           edgecolor="black", linewidth=0.3, zorder=3)
            ax.axvline(EXPECTED[comp][sp], ls="--", lw=1.1, color=COLORS[sp], zorder=1)
        ax.set_title(comp, fontsize=11)
        ax.set_xlim(*x_lim)
        ax.grid(True, axis="x", linewidth=0.3)
    axes[-1].set_xlabel("Ratio (log scale) — dashed = design-expected")
    fig.suptitle(title, fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(out_png, dpi=200)
    fig.savefig(out_png.replace(".png", ".pdf"))
    plt.close(fig)
    print(f"  wrote {out_png}")


if __name__ == "__main__":
    rep = load_reported()
    rec = load_recalculated()
    rep.to_csv(os.path.join(RESULTS, "reported_ratios_long.csv"), index=False)
    rec.to_csv(os.path.join(RESULTS, "recalculated_ratios_long.csv"), index=False)
    print(f"reported: {len(rep)} points, {rep['label'].nunique()} entries -> {sorted(rep['label'].unique())}")
    print(f"recalc  : {len(rec)} points, {rec['label'].nunique()} entries -> {sorted(rec['label'].unique())}")
    boxwhisker(rep, "Participant-reported ratios (n=%d entries)" % rep["label"].nunique(),
               os.path.join(FIGURES, "FIG04A_reported_ratios.png"))
    boxwhisker(rec, "Ratios recalculated from peptide quant (n=%d entries)" % rec["label"].nunique(),
               os.path.join(FIGURES, "FIG04B_recalculated_ratios.png"))
