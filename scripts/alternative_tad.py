import argparse
import subprocess
from pathlib import Path
import cooler
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

from analyze_opt import ROOT, DATA, OUT, SAMPLES, RES, CHROM, START, END, WINDOW, savefig, overlap_region

CONTEXT_START, CONTEXT_END = 7_500_000, 13_500_000


def match_boundaries(a, b, tolerance=RES):

    a, b = sorted(set(a)), sorted(set(b))
    i = j = 0
    pairs = []
    while i < len(a) and j < len(b):
        if abs(a[i] - b[j]) <= tolerance:
            pairs.append((a[i], b[j], abs(a[i] - b[j])))
            i += 1
            j += 1
        elif a[i] < b[j]:
            i += 1
        else:
            j += 1
    return pairs


def boundary_positions(domains):


    return sorted({int(x) for x in list(domains.start) + list(domains.end - RES) if START <= x < END})


def best_iou(left, right, other):
    if other.empty:
        return 0.0
    intersection = np.maximum(0, np.minimum(right, other.end) - np.maximum(left, other.start))
    union = (right - left) + (other.end - other.start) - intersection
    return float((intersection / union).max())


def main(binary):
    comparison, details, domain_comparison = [], [], []
    for name, accession in SAMPLES.items():
        c = cooler.Cooler(f"{DATA / (accession + '.mcool')}::/resolutions/{RES}")
        matrix = c.matrix(balance=False).fetch((CHROM, CONTEXT_START, CONTEXT_END))
        assert matrix.shape == ((CONTEXT_END-CONTEXT_START)//RES,) * 2
        assert np.array_equal(matrix, matrix.T) and (matrix >= 0).all()
        matrix_path = OUT / f"{name}_ontad_input.tsv"
        np.savetxt(matrix_path, matrix, fmt="%d", delimiter="\t")
        prefix = OUT / f"{name}_ontad"
        command = [str(binary), str(matrix_path), "-penalty", "0.1", "-minsz", "3", "-maxsz", "80", "-ldiff", "1.96", "-lsize", "5", "-o", str(prefix)]
        subprocess.run(command, stdout=subprocess.DEVNULL, check=True)
        matrix_path.unlink()
        calls = pd.read_csv(str(prefix) + ".tad", sep=r"\s+", header=None,
                            names=["start_bin", "end_bin", "level", "mean", "score"])
        calls = calls[calls.level > 0].copy()
        if not calls.empty:
            assert (calls.start_bin >= 1).all() and (calls.end_bin <= len(matrix)).all()

        calls["chrom"] = CHROM
        calls["start"] = CONTEXT_START + (calls.start_bin - 1) * RES
        calls["end"] = CONTEXT_START + calls.end_bin * RES
        calls = calls[["chrom", "start", "end", "level", "mean", "score", "start_bin", "end_bin"]]
        calls.to_csv(OUT / f"{name}_ontad_domains.tsv", sep="\t", index=False)
        local = overlap_region(calls)
        ins = pd.read_csv(OUT / f"{name}_insulation_locus.tsv", sep="\t")
        a = ins.loc[ins[f"is_boundary_{WINDOW}"], "start"].to_numpy(int)
        for subset, domains in [("level1", local[local.level == 1]), ("all_levels", local)]:
            endpoints = boundary_positions(domains)
            pairs = match_boundaries(a, endpoints)
            for x, y, d in pairs:
                details.append({"sample": name, "ontad_subset": subset, "cooltools_bp": x, "ontad_bp": y, "distance_bp": d})
            recall = len(pairs)/len(a) if len(a) else np.nan
            precision = len(pairs)/len(endpoints) if endpoints else np.nan
            comparison.append({"sample": name, "ontad_subset": subset, "cooltools_boundaries": len(a),
                "ontad_boundaries": len(endpoints), "ontad_domains_overlapping": len(domains),
                "matched_one_to_one": len(pairs), "tolerance_bp": RES,
                "fraction_cooltools_matched": recall, "fraction_ontad_matched": precision,
                "boundary_jaccard": len(pairs)/(len(a)+len(endpoints)-len(pairs)) if len(a)+len(endpoints) else np.nan})
            reference = overlap_region(pd.read_csv(OUT / f"{name}_tads_{WINDOW}.tsv", sep="\t"))
            for d in reference.itertuples():
                domain_comparison.append({"sample": name, "ontad_subset": subset, "start": d.start, "end": d.end,
                                          "best_domain_iou": best_iou(d.start, d.end, domains)})
        endpoints = boundary_positions(local)
        fig, (ax, profile) = plt.subplots(2, 1, figsize=(8, 9), gridspec_kw={"height_ratios": [3, 1]})
        crop = c.matrix(balance=False).fetch((CHROM, START, END))
        ax.imshow(np.log1p(crop), extent=[START/1e6, END/1e6, START/1e6, END/1e6], origin="lower", cmap="magma_r")
        for d in local.itertuples():
            ax.add_patch(Rectangle((d.start/1e6, d.start/1e6), (d.end-d.start)/1e6, (d.end-d.start)/1e6,
                fill=False, edgecolor="#1b9e77", lw=1.5 if d.level == 1 else 0.7))
        for x in a:
            ax.axvline(x/1e6, color="#2166ac", ls="--", alpha=0.6)
            ax.axhline(x/1e6, color="#2166ac", ls="--", alpha=0.6)
        ax.set(xlim=(START/1e6, END/1e6), ylim=(START/1e6, END/1e6), xlabel="chr17 (Mb)", ylabel="chr17 (Mb)",
               title=f"{name}: OnTAD domains (green); cooltools boundaries (blue)")
        profile.plot(ins.start/1e6, ins[f"log2_insulation_score_{WINDOW}"], color="black")
        for x in a:
            profile.axvline(x/1e6, color="#2166ac", ls="--")
        for x in endpoints:
            profile.axvline(x/1e6, color="#1b9e77", ls=":", alpha=0.7)
        profile.set(xlim=(START/1e6, END/1e6), xlabel="chr17 (Mb)", ylabel="log2 insulation")
        savefig(fig, f"{name}_ontad_comparison.png")
    pd.DataFrame(comparison).to_csv(OUT / "alternative_tad_comparison.tsv", sep="\t", index=False)
    pd.DataFrame(details, columns=["sample", "ontad_subset", "cooltools_bp", "ontad_bp", "distance_bp"]).to_csv(OUT / "boundary_matches.tsv", sep="\t", index=False)
    pd.DataFrame(domain_comparison, columns=["sample", "ontad_subset", "start", "end", "best_domain_iou"]).to_csv(OUT / "domain_comparison.tsv", sep="\t", index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--ontad", type=Path, default=ROOT / "tools" / "OnTAD" / "src" / "OnTAD")
    main(parser.parse_args().ontad.resolve())
