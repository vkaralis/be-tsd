"""Command-line simulation with JSON output and full scenario metadata."""
import argparse
import json
from pathlib import Path

from . import METHODS, Design, simulate


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--method", choices=METHODS, required=True)
    parser.add_argument("--n1", type=int, default=24)
    parser.add_argument("--gmr", type=float, default=0.95)
    parser.add_argument("--cv", type=float, default=0.2, help="CV fraction, e.g. 0.20")
    parser.add_argument("--upper-limit", type=int)
    parser.add_argument("--assumed-gmr", type=float)
    parser.add_argument("--power-engine", choices=["potvin", "legacy"])
    parser.add_argument("--trials", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=2026)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    try:
        design = Design(args.method, args.upper_limit, args.assumed_gmr,
                        power_engine=args.power_engine)
        result = simulate(design, args.n1, args.gmr, args.cv, args.trials, args.seed)
    except (ValueError, RuntimeError) as error:
        parser.error(str(error))
    payload = json.dumps(result, indent=2, allow_nan=False) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(payload, encoding="utf-8")
    print(payload, end="")


if __name__ == "__main__":
    main()
