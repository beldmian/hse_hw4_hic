from pathlib import Path
import json
import subprocess
import sys
import importlib.metadata

import cooler
import cooltools
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

ROOT = Path(__file__).resolve().parents[1]
DATA, OUT = ROOT / "data", ROOT / "results"
SAMPLES = {"MCG037": "4DNFIJ2JKO7D", "MCG011": "4DNFI8I9LN74"}
RES = 25_000
CHROM, START, END = "chr17", 10_000_000, 11_000_000
REGION = f"{CHROM}:{START}-{END}"
WINDOWS = [100_000, 200_000, 300_000, 500_000, 1_000_000]
WINDOW = 200_000
COLORS = {"MCG037": "#2166ac", "MCG011": "#b2182b"}


def savefig(fig, name):
    fig.tight_layout()
    fig.savefig(OUT / name, dpi=180)
    plt.close(fig)


def write_boundary_bed(boundaries, path, strength_column):

    bed = boundaries[["chrom", "start", "end"]].copy()
    bed["name"] = [f"boundary_{i+1}" for i in range(len(bed))]
    strength = boundaries[strength_column].to_numpy(float)
    if not np.isfinite(strength).all() or np.any(strength < 0):
        raise ValueError("Boundary strength must be finite and nonnegative")

    bed["score"] = np.rint(1000 * strength / (1 + strength)).astype(int)
    bed["strand"] = "."
    bed.to_csv(path, sep="\t", header=False, index=False)


def contact_distance(c, view):

    kwargs = dict(view_df=view, ignore_diags=2, chunksize=1_000_000, nproc=1)
    balanced = cooltools.expected_cis(c, clr_weight_name="weight", **kwargs)
    raw = cooltools.expected_cis(c, clr_weight_name=None, **kwargs)
    keys = ["region1", "region2", "dist"]


    return balanced.drop(columns=["count.sum", "count.avg"]).merge(
        raw[keys + ["count.sum", "count.avg"]], on=keys, validate="one_to_one")


def domains_from_boundaries(table, window):


    boundaries = table.loc[table[f"is_boundary_{window}"].fillna(False), "start"].to_numpy(int)
    bad = table["is_bad_bin"].to_numpy(bool) | ~np.isfinite(table[f"log2_insulation_score_{window}"])
    rows = []
    starts = table.start.to_numpy()
    for left, right in zip(boundaries[:-1], boundaries[1:]):
        if not bad[(starts >= left) & (starts <= right)].any():
            rows.append((CHROM, left, right))
    return pd.DataFrame(rows, columns=["chrom", "start", "end"])


def overlap_region(table):
    return table[(table.start < END) & (table.end > START)].copy()


def main():
    OUT.mkdir(exist_ok=True)
    summary = {"parameters": {"resolution_bp": RES, "region": REGION, "chromosome": CHROM,
                "windows_bp": WINDOWS, "main_window_bp": WINDOW, "insulation_view": "whole chr17",
                "threshold": "Li", "ignore_diags": 2, "min_frac_valid_pixels": 0.66},
               "versions": {p: importlib.metadata.version(p) for p in ["cooler", "cooltools", "numpy", "pandas", "scipy", "matplotlib"]},
               "samples": {}}
    curves, profiles, matrices, counts = {}, {}, {}, []
    for name, accession in SAMPLES.items():
        path = DATA / f"{accession}.mcool"
        uri = f"{path}::/resolutions/{RES}"
        c = cooler.Cooler(uri)
        if "weight" not in c.bins().columns:
            raise ValueError(f"{accession}: missing balancing weights")
        info = dict(c.info)
        view = pd.DataFrame([[CHROM, 0, int(c.chromsizes[CHROM]), CHROM]], columns=["chrom", "start", "end", "name"])
        bins = c.bins().fetch(REGION)
        bins.to_csv(OUT / f"{name}_bins.tsv", sep="\t", index=False)

        matrix = c.matrix(balance="weight").fetch(REGION)
        matrices[name] = matrix
        pixels = c.matrix(balance="weight", as_pixels=True, join=True).fetch(REGION)
        pixels.to_csv(OUT / f"{name}_contacts.tsv", sep="\t", index=False)
        cmd = [sys.executable, "-m", "cooler", "dump", "--header", "--join", "--balanced", "--range", REGION, uri]
        with (OUT / f"{name}_cooler_dump_balanced.tsv").open("w") as out:
            subprocess.run(cmd, stdout=out, check=True)
        print(f"{name}: expected cis contacts", flush=True)
        expected = contact_distance(c, view)
        expected.to_csv(OUT / f"{name}_contacts_vs_distance.tsv", sep="\t", index=False)
        curves[name] = expected
        print(f"{name}: insulation on entire chr17, all five windows", flush=True)
        ins = cooltools.insulation(c, WINDOWS, view_df=view, ignore_diags=2, threshold="Li",
                min_frac_valid_pixels=0.66, append_raw_scores=True, chunksize=1_000_000)
        ins.to_csv(OUT / f"{name}_insulation_chr17.tsv", sep="\t", index=False)
        locus = overlap_region(ins)
        locus.to_csv(OUT / f"{name}_insulation_locus.tsv", sep="\t", index=False)
        profiles[name] = locus
        for w in WINDOWS:
            domains = domains_from_boundaries(ins, w)
            domains.to_csv(OUT / f"{name}_tads_{w}.tsv", sep="\t", index=False)
            local_domains = overlap_region(domains)
            fully_inside = domains[(domains.start >= START) & (domains.end <= END)]
            nboundary = int(locus[f"is_boundary_{w}"].fillna(False).sum())
            counts.append({"sample": name, "window_bp": w, "boundaries_in_locus": nboundary,
                "tads_overlapping_locus": len(local_domains), "tads_fully_inside_locus": len(fully_inside),
                "tads_chr17": len(domains), "finite_score_bins_locus": int(np.isfinite(locus[f"log2_insulation_score_{w}"]).sum())})
            if w == WINDOW:
                boundaries = locus[locus[f"is_boundary_{w}"].fillna(False)]
                write_boundary_bed(boundaries, OUT / f"{name}_boundaries.bed", f"boundary_strength_{w}")
        row100 = expected[expected.dist_bp == 100_000].iloc[0]
        summary["samples"][name] = {"accession": accession, "info": info,
                "resolutions": cooler.fileops.list_coolers(path), "chromnames": list(c.chromnames),
                "bins_columns": list(c.bins().columns), "pixel_columns": list(c.pixels().columns),
                "contact_table_columns": list(pixels.columns), "locus_bins": len(bins),
                "valid_weight_bins": int(np.isfinite(bins.weight).sum()),
                "raw_mean_at_100kb": float(row100["count.avg"]), "balanced_mean_at_100kb": float(row100["balanced.avg"])}
    pd.DataFrame(counts).to_csv(OUT / "window_counts.tsv", sep="\t", index=False)
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2, default=int))

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.8), sharex=True, sharey=True, layout="constrained")
    values = np.concatenate([m[np.isfinite(m) & (m > 0)] for m in matrices.values()])
    norm = LogNorm(vmin=np.quantile(values, 0.05), vmax=np.quantile(values, 0.98))
    for ax, (name, m) in zip(axes, matrices.items()):
        im = ax.imshow(np.ma.masked_invalid(np.where(m > 0, m, np.nan)), norm=norm, cmap="magma_r",
                       origin="lower", extent=[START/1e6, END/1e6, START/1e6, END/1e6])
        ax.set(title=f"{name} / {SAMPLES[name]}", xlabel="chr17 (Mb)", ylabel="chr17 (Mb)")
    fig.colorbar(im, ax=axes, label="Balanced contact (shared log scale)", shrink=0.75)
    fig.savefig(OUT / "contact_maps.png", dpi=180)
    plt.close(fig)

    fig, axes = plt.subplots(1, 3, figsize=(14, 4))
    for name, e in curves.items():
        for ax, col in zip(axes[:2], ["count.avg", "balanced.avg"]):
            z = e[(e.dist_bp > 0) & np.isfinite(e[col]) & (e[col] > 0)]
            ax.loglog(z.dist_bp, z[col], color=COLORS[name], label=name)
        z = e[(e.dist_bp > 0) & np.isfinite(e["balanced.avg"]) & (e["balanced.avg"] > 0)]
        axes[2].loglog(z.dist_bp, z["balanced.avg"] / summary["samples"][name]["balanced_mean_at_100kb"], label=name, color=COLORS[name])
    for ax, title in zip(axes, ["Raw mean contacts", "Balanced mean contacts", "Balanced / value at 100 kb"]):
        ax.set(xlabel="Genomic separation (bp)", ylabel=title, title=title)
        ax.grid(alpha=0.2, which="both")
        ax.legend()
    savefig(fig, "contacts_vs_distance_comparison.png")

    fig, axes = plt.subplots(2, 1, figsize=(10, 6), sharex=True, sharey=True)
    for ax, (name, locus) in zip(axes, profiles.items()):
        score = f"log2_insulation_score_{WINDOW}"
        b = locus[locus[f"is_boundary_{WINDOW}"].fillna(False)]
        ax.plot(locus.start / 1e6, locus[score], label=name, color=COLORS[name])
        for x in b.start:
            ax.axvline(x / 1e6, color=COLORS[name], linestyle="--", alpha=0.7)
        ax.scatter(b.start / 1e6, b[score], color=COLORS[name], zorder=3)
        ax.set(ylabel="log2 insulation", title=f"{name}: {WINDOW//1000} kb window, dashed lines = called boundaries")
        ax.grid(alpha=0.2)
        ax.legend()
    axes[-1].set_xlabel("chr17 (Mb)")
    savefig(fig, "insulation_comparison.png")

    fig, axes = plt.subplots(1, 3, figsize=(13, 4))
    df = pd.DataFrame(counts)
    for name in SAMPLES:
        sub = df[df["sample"] == name]
        for ax, col in zip(axes, ["boundaries_in_locus", "tads_overlapping_locus", "tads_fully_inside_locus"]):
            ax.plot(sub.window_bp / 1000, sub[col], "o-", label=name, color=COLORS[name])
            ax.set(xlabel="Window (kb)", ylabel="Count", title=col.replace("_", " "))
            ax.grid(alpha=0.2)
            ax.yaxis.get_major_locator().set_params(integer=True)
    axes[0].legend()
    savefig(fig, "window_vs_tad_count.png")
    print("Analysis complete", flush=True)


if __name__ == "__main__":
    main()
