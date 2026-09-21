import hashlib
import json
import sys
import tempfile
from pathlib import Path
import cooler
import numpy as np
import pandas as pd
from analyze_opt import DATA, OUT, SAMPLES, RES, CHROM, START, END, WINDOWS, WINDOW, domains_from_boundaries, overlap_region, contact_distance
from alternative_tad import match_boundaries, best_iou, boundary_positions


def main():
    assert match_boundaries([100, 110], [105], 10) == [(100, 105, 5)]
    assert match_boundaries([], [105], 10) == []
    assert len(match_boundaries([100, 200], [125, 225], 25)) == 2
    assert match_boundaries([100], [126], 25) == []
    synthetic = pd.DataFrame({"start": [0, 25, 50, 75], "is_bad_bin": [False]*4,
        "log2_insulation_score_100": [0.0]*4, "is_boundary_100": [True, False, True, True]})
    assert len(domains_from_boundaries(synthetic, 100)) == 2
    synthetic.loc[1, "is_bad_bin"] = True
    assert len(domains_from_boundaries(synthetic, 100)) == 1
    synthetic.loc[3, "log2_insulation_score_100"] = np.nan
    assert domains_from_boundaries(synthetic, 100).empty
    assert best_iou(0, 100, pd.DataFrame({"start": [50], "end": [150]})) == 1/3
    shared = pd.DataFrame({"start": [10_000_000, 10_500_000], "end": [10_525_000, 11_025_000]})
    assert boundary_positions(shared) == [10_000_000, 10_500_000]
    with tempfile.TemporaryDirectory() as tmp:
        path = str(Path(tmp) / "missing_weight.cool")
        bins = pd.DataFrame({"chrom": [CHROM]*12, "start": np.arange(12)*RES, "end": np.arange(1,13)*RES, "weight": 1.0})
        bins.loc[5, "weight"] = np.nan
        i, j = np.triu_indices(12)
        pixels = pd.DataFrame({"bin1_id": i, "bin2_id": j, "count": np.where((i == 5) | (j == 5), 1000, 10)})
        cooler.create_cooler(path, bins, pixels)
        view = pd.DataFrame([[CHROM, 0, 12*RES, CHROM]], columns=["chrom", "start", "end", "name"])
        result = contact_distance(cooler.Cooler(path), view)
        row = result[result.dist == 4].iloc[0]
        assert np.isclose(row["count.avg"], pixels.loc[j-i == 4, "count"].sum()/8)
        assert np.isclose(row["balanced.avg"], 10.0)
    if "--algorithms-only" in sys.argv:
        print("Algorithm checks passed")
        return
    counts = pd.read_csv(OUT / "window_counts.tsv", sep="\t")
    summary = json.loads((OUT / "summary.json").read_text())
    assert len(counts) == 2 * len(WINDOWS)
    for name, accession in SAMPLES.items():
        metadata = json.loads((DATA / "sources.json").read_text())[accession]
        path = DATA / f"{accession}.mcool"
        assert path.stat().st_size == metadata["file_size"]
        with path.open("rb") as f:
            assert hashlib.file_digest(f, "md5").hexdigest() == metadata["md5sum"]
        info = summary["samples"][name]
        assert info["locus_bins"] == (END-START)//RES
        assert info["valid_weight_bins"] == info["locus_bins"], "Choose a well-covered locus; do not hide missing weights"
        c = cooler.Cooler(f"{path}::/resolutions/{RES}")
        matrix = c.matrix(balance="weight").fetch((CHROM, START, END))
        assert matrix.shape == ((END-START)//RES,) * 2
        assert np.allclose(matrix, matrix.T, equal_nan=True)
        pixels = pd.read_csv(OUT / f"{name}_contacts.tsv", sep="\t")
        bins = pd.read_csv(OUT / f"{name}_bins.tsv", sep="\t").set_index("start")
        expected = pixels["count"] * pixels.start1.map(bins.weight) * pixels.start2.map(bins.weight)
        assert np.allclose(pixels.balanced, expected, equal_nan=True)
        balanced = pd.read_csv(OUT / f"{name}_cooler_dump_balanced.tsv", sep="\t")
        assert len(balanced) == len(pixels) > 0
        assert np.array_equal(balanced["count"], pixels["count"])
        assert np.allclose(balanced.balanced, pixels.balanced, rtol=1e-5, equal_nan=True)
        table = pd.read_csv(OUT / f"{name}_insulation_chr17.tsv", sep="\t")
        locus = overlap_region(table)
        for w in WINDOWS:
            saved = pd.read_csv(OUT / f"{name}_tads_{w}.tsv", sep="\t")
            generated = domains_from_boundaries(table, w)
            pd.testing.assert_frame_equal(saved, generated, check_dtype=False)
            row = counts[(counts["sample"] == name) & (counts.window_bp == w)].iloc[0]
            assert row.tads_overlapping_locus == len(overlap_region(generated))
            assert row.tads_fully_inside_locus == len(generated[(generated.start >= START) & (generated.end <= END)])
            assert row.boundaries_in_locus == int(locus[f"is_boundary_{w}"].sum())
            assert row.finite_score_bins_locus == len(locus)
        boundary_rows = locus[locus[f"is_boundary_{WINDOW}"]]
        bed_path = OUT / f"{name}_boundaries.bed"
        if len(boundary_rows):
            bed = pd.read_csv(bed_path, sep="\t", header=None)
            assert len(bed) == len(boundary_rows) and bed.shape[1] == 6
            assert (bed[0] == CHROM).all() and (bed[2]-bed[1] == RES).all()
            assert ((bed[4] >= 0) & (bed[4] <= 1000)).all()
            strength = boundary_rows[f"boundary_strength_{WINDOW}"].to_numpy()
            assert np.array_equal(bed[4], np.rint(1000*strength/(1+strength)).astype(int))
        else:
            assert bed_path.stat().st_size == 0
        calls = pd.read_csv(OUT / f"{name}_ontad_domains.tsv", sep="\t")
        assert (calls.level > 0).all() and (calls.end > calls.start).all()
        assert (calls.start == 7_500_000 + (calls.start_bin-1)*RES).all()
        assert (calls.end == 7_500_000 + calls.end_bin*RES).all()
        assert (OUT.parent / "report" / "higlass" / f"{name}_HiGlass.png").stat().st_size > 10000
    comparison = pd.read_csv(OUT / "alternative_tad_comparison.tsv", sep="\t")
    assert (comparison.matched_one_to_one <= comparison.cooltools_boundaries).all()
    assert (comparison.matched_one_to_one <= comparison.ontad_boundaries).all()
    print("PASS: source MD5, balanced contacts, CLI dumps, five windows, TAD counts, BED strengths, OnTAD, HiGlass files")


if __name__ == "__main__":
    main()
