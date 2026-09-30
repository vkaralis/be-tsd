"""Plot acceptance and mean total sample size from run_study.py CSV output."""
import argparse
import csv
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("--output-dir", type=Path, default=Path("results/plots"))
    args = parser.parse_args()
    groups = defaultdict(list)
    with args.input.open(newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream):
            groups[(row["method"], row["true_cv"], row["n1"])].append(row)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for (method, cv, n1), rows in groups.items():
        fig, axes = plt.subplots(1, 2, figsize=(10, 4), layout="constrained")
        for ul in sorted({r["upper_limit"] for r in rows}):
            curve = sorted((r for r in rows if r["upper_limit"] == ul), key=lambda r: float(r["true_gmr"]))
            x = [float(r["true_gmr"]) for r in curve]
            label = f"UL={ul}" if ul else "No UL"
            y = [100*float(r["acceptance_probability"]) for r in curve]
            line, = axes[0].plot(x, y, "o-", label=label)
            axes[0].fill_between(x, [100*float(r["mc_ci_low"]) for r in curve],
                                 [100*float(r["mc_ci_high"]) for r in curve],
                                 alpha=.15, color=line.get_color())
            axes[1].plot(x, [float(r["mean_total_n"]) for r in curve], "o-", label=label)
        axes[0].set(ylabel="BE acceptance (%)", ylim=(0, 100))
        axes[1].set(ylabel="Mean total sample size")
        for ax in axes:
            ax.set_xlabel("True GMR (Test / Reference)")
            ax.grid(alpha=.25)
            ax.legend()
        fig.suptitle(f"{method} | CVw={float(cv):.0%} | N1={n1}")
        fig.savefig(args.output_dir / f"{method}_cv{cv}_n{n1}.png", dpi=180)
        plt.close(fig)
    # Figure 4 compares curves with different N1, so plot it separately.
    figure4_curves = []
    for n1, ul in [(72, 72), (48, 150)]:
        curve = [r for r in groups.get(("TSD-2", "0.2", str(n1)), [])
                 if r["upper_limit"] == str(ul)]
        figure4_curves.append((n1, ul, sorted(curve, key=lambda r: float(r["true_gmr"]))))
    if all(curve for _, _, curve in figure4_curves):
        fig, ax = plt.subplots(figsize=(6, 4), layout="constrained")
        for n1, ul, curve in figure4_curves:
            ax.plot([float(r["true_gmr"]) for r in curve],
                    [100*float(r["acceptance_probability"]) for r in curve],
                    "o-", label=f"N1={n1}, UL={ul}")
        ax.set(xlabel="True GMR (Test / Reference)", ylabel="BE acceptance (%)",
               ylim=(0, 100), title="Figure 4 scenarios | TSD-2 | CVw=20%")
        ax.grid(alpha=.25)
        ax.legend()
        fig.savefig(args.output_dir / "figure4_comparison.png", dpi=180)
        plt.close(fig)


if __name__ == "__main__":
    main()
