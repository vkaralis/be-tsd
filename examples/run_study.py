"""Run article grids or a comparison of the five methods; save CSV.

Usage (after installing the package):
    python examples/run_study.py --preset comparison --trials 1000
    python examples/run_study.py --preset tables --trials 1000000
    python examples/run_study.py --preset figures --trials 100000
"""
import argparse
import csv
from itertools import product
import json
from pathlib import Path

import numpy as np

from twostage_be import Design, simulate


def scenarios(preset):
    if preset == "tables":
        for method, cv, n1, ul in product(["TSD-1", "TSD-2"], [.2, .4, .6],
                                         [12, 24, 36, 60], [62, 100, 150, 500]):
            yield Design(method, ul), n1, 1.25, cv
    elif preset == "figures":
        gmrs = [round(1 + i * .025, 3) for i in range(11)]
        for method, cv, n1, ul, gmr in product(["TSD-1", "TSD-2"], [.2, .4],
                                              [24, 48], [62, 150], gmrs):
            yield Design(method, ul), n1, gmr, cv
        # Figure 4: the N1=48/UL=150 curve is already in the above grid.
        for gmr in gmrs:
            yield Design("TSD-2", 72), 72, gmr, .2
    else:
        for method, cv, gmr in product(["TSD-1", "TSD-2", "B", "C", "D"],
                                       [.2, .4], [.95, 1.0, 1.25]):
            # Finite UL is mandatory for TSD and optional for Potvin.
            yield Design(method, 150 if method.startswith("TSD") else None), 24, gmr, cv


def references():
    path = Path(__file__).resolve().parents[1] / "docs" / "article_tables.csv"
    with path.open(newline="", encoding="utf-8") as stream:
        return {(r["method"], float(r["true_cv"]), int(r["n1"]), int(r["upper_limit"])):
                float(r["published_percent"]) for r in csv.DictReader(stream)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--preset", choices=["comparison", "tables", "figures"], default="comparison")
    parser.add_argument("--trials", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=2026)
    parser.add_argument("--output", type=Path, default=Path("results/study.csv"))
    args = parser.parse_args()
    if args.trials < 1:
        parser.error("trials must be positive")
    settings = list(scenarios(args.preset))
    seeds = np.random.SeedSequence(args.seed).spawn(len(settings))
    refs = references()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as stream:
        writer = None
        for i, ((design, n1, gmr, cv), child) in enumerate(zip(settings, seeds), 1):
            seed = int(child.generate_state(1)[0])
            result = simulate(design, n1, gmr, cv, args.trials, seed)
            result["reasons"] = json.dumps(result["reasons"], sort_keys=True)
            published = refs.get((design.method, cv, n1, design.upper_limit)) if args.preset == "tables" else None
            result["published_percent"] = published
            result["difference_percentage_points"] = (100 * result["acceptance_probability"] - published
                                                      if published is not None else None)
            if writer is None:
                writer = csv.DictWriter(stream, fieldnames=list(result))
                writer.writeheader()
            writer.writerow(result)
            stream.flush()
            print(f"{i}/{len(settings)} {design.method} N1={n1} CV={cv} GMR={gmr} UL={design.upper_limit}",
                  flush=True)
    args.output.with_suffix(".metadata.json").write_text(json.dumps({
        "preset": args.preset, "master_seed": args.seed, "trials_per_scenario": args.trials,
        "scenario_count": len(settings), "seed_strategy": "SeedSequence.spawn in scenario order",
    }, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
