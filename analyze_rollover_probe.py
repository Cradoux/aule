#!/usr/bin/env python3
import argparse
import json
import os
from typing import Dict, Any

import pandas as pd
import numpy as np

try:
    import matplotlib.pyplot as plt
except Exception:
    plt = None


def basic_tallies(df: pd.DataFrame) -> Dict[str, Any]:
    hopped = (df.f1 != df.f0)
    out = {
        "rows": int(len(df)),
        "hops": int(hopped.sum()),
        "hop_rate": float(hopped.mean()),
        "kneg_counts": df.kneg.value_counts().sort_index().to_dict(),
    }
    # Validate kneg domain
    uniq = set(df.kneg.unique().tolist())
    invalid = sorted([int(v) for v in uniq if v not in {0, 1, 2}])
    if invalid:
        out["invalid_kneg_values"] = invalid
    # Dimensions (if full frame present)
    try:
        w = int(df.x.max()) + 1
        h = int(df.y.max()) + 1
        out["width"] = w
        out["height"] = h
    except Exception:
        pass
    return out


def invariant_checks(df: pd.DataFrame) -> Dict[str, Any]:
    bad_gpu = (df.face_gpu != df.f1)
    bad_nf = (df.f1 != df.nf) & (df.f1 != df.f0)
    return {
        "gpu_vs_f1_mismatch": int(bad_gpu.sum()),
        "illegal_final_face": int(bad_nf.sum()),
    }


def hop_predicate_checks(df: pd.DataFrame, eps: float = 1e-6) -> Dict[str, Any]:
    if not set(["wa", "wb", "wc"]).issubset(df.columns):
        return {"note": "wa/wb/wc columns not present; skipping predicate checks"}
    wmin = df[["wa", "wb", "wc"]].min(axis=1)
    neg = wmin < -eps
    hop = (df.f1 != df.f0)
    return {
        "hop_but_not_negative": int((hop & (~neg)).sum()),
        "no_hop_but_negative": int(((~hop) & neg).sum()),
    }


def neighbor_mapping_sample(df: pd.DataFrame, limit: int = 20) -> pd.DataFrame:
    wrong_neighbor = (df.f1 != df.f0) & (df.f1 != df.nf)
    cols = ["x", "y", "f0", "kneg", "nf", "f1"]
    cols = [c for c in cols if c in df.columns]
    return df.loc[wrong_neighbor, cols].head(limit).copy()


def plot_hops(df: pd.DataFrame, out_path: str) -> bool:
    if plt is None:
        return False
    if not set(["x", "y", "f0", "f1"]).issubset(df.columns):
        return False
    try:
        w = int(df.x.max()) + 1
        h = int(df.y.max()) + 1
        hop_mask = np.zeros((h, w), dtype=np.uint8)
        hop_df = df.loc[df.f1 != df.f0, ["x", "y"]]
        hop_mask[hop_df.y.to_numpy(), hop_df.x.to_numpy()] = 1
        plt.figure(figsize=(10, 5), dpi=150)
        plt.imshow(hop_mask, cmap="magma", origin="upper", interpolation="nearest")
        plt.colorbar(label="hop (1=yes)")
        plt.title("Rollover hops (f1 != f0)")
        plt.xlabel("x")
        plt.ylabel("y")
        plt.tight_layout()
        plt.savefig(out_path)
        plt.close()
        return True
    except Exception:
        return False


def main() -> None:
    ap = argparse.ArgumentParser(description="Analyze rollover_probe.csv seams and invariants")
    ap.add_argument("--csv", default="rollover_probe.csv", help="Path to rollover_probe.csv")
    ap.add_argument("--eps", type=float, default=1e-6, help="Rollover epsilon for negativity test")
    ap.add_argument("--out_json", default="rollover_probe_summary.json", help="Summary JSON output path")
    ap.add_argument("--out_img", default="rollover_probe_hops.png", help="Heat image output path")
    args = ap.parse_args()

    if not os.path.exists(args.csv):
        raise SystemExit(f"CSV not found: {args.csv}")

    # Read with minimal dtype fuss; let pandas infer
    df = pd.read_csv(args.csv)

    required = ["x", "y", "f0", "kneg", "nf", "f1", "face_gpu", "tri_gpu"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        print({"warning": f"Missing columns: {missing}"})

    results: Dict[str, Any] = {}
    results["tallies"] = basic_tallies(df)
    results["invariants"] = invariant_checks(df)
    results["predicate"] = hop_predicate_checks(df, eps=args.eps)

    # Sample wrong neighbor rows
    sample = neighbor_mapping_sample(df, limit=20)
    if len(sample) > 0:
        results["wrong_neighbor_sample_rows"] = int(len(sample))
        sample.to_csv("rollover_probe_wrong_neighbor_sample.csv", index=False)
        results["wrong_neighbor_sample_csv"] = "rollover_probe_wrong_neighbor_sample.csv"

    # Plot hops heat
    results["heat_saved"] = plot_hops(df, args.out_img)
    results["heat_path"] = args.out_img if results["heat_saved"] else None

    # Print and write JSON summary
    print(json.dumps(results, indent=2))
    with open(args.out_json, "w", encoding="utf-8") as f:
        json.dump(results, f, indent=2)


if __name__ == "__main__":
    main()



