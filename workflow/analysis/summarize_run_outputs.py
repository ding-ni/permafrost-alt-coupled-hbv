# -*- coding: utf-8 -*-
import argparse
import json
from pathlib import Path


def load_json(path):
    with Path(path).open("r", encoding="utf-8") as handle:
        return json.load(handle)


def main():
    parser = argparse.ArgumentParser(description="Summarize standardized HBV output files.")
    parser.add_argument("--original-summary", default="", help="Path to original_hbv_summary.json")
    parser.add_argument("--improved-metadata", default="", help="Path to improved HBV metadata.json")
    parser.add_argument("--alt-summary", default="", help="Path to ALT-coupled summary.json")
    args = parser.parse_args()

    if args.original_summary:
        summary = load_json(args.original_summary)
        print("[original_hbv]")
        print(json.dumps(summary, ensure_ascii=False, indent=2))

    if args.improved_metadata:
        metadata = load_json(args.improved_metadata)
        reduced = {
            "run_dir": metadata.get("run_dir"),
            "optimized_params": metadata.get("optimized_params"),
            "metrics_objective_window": metadata.get("metrics_objective_window"),
            "metrics_full_period": metadata.get("metrics_full_period"),
        }
        print("[improved_hbv]")
        print(json.dumps(reduced, ensure_ascii=False, indent=2))

    if args.alt_summary:
        summary = load_json(args.alt_summary)
        print("[alt_coupled_hbv]")
        print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
