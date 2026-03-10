# -*- coding: utf-8 -*-
import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


def main():
    parser = argparse.ArgumentParser(description="Plot a simulated hydrograph against observations.")
    parser.add_argument("--timeseries", required=True, help="Path to a *_timeseries.csv file.")
    parser.add_argument("--output", default="", help="Optional output PNG path.")
    args = parser.parse_args()

    df = pd.read_csv(args.timeseries, parse_dates=["date"])
    fig, ax = plt.subplots(figsize=(10, 4))
    if "q_obs" in df:
        ax.plot(df["date"], df["q_obs"], label="Observed", linewidth=1.4)
    if "q_sim" in df:
        ax.plot(df["date"], df["q_sim"], label="Simulated", linewidth=1.2)
    ax.set_xlabel("Date")
    ax.set_ylabel("Discharge")
    ax.legend()
    ax.set_title(Path(args.timeseries).stem)
    fig.tight_layout()

    if args.output:
        fig.savefig(args.output, dpi=150)
    else:
        plt.show()


if __name__ == "__main__":
    main()
