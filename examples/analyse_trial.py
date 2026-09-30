"""Analyse supplied log differences and follow the interim decision."""
import numpy as np
from twostage_be import Design, analyse, interim


def main():
    design = Design("TSD-2", upper_limit=150)
    # Replace with your own complete data: one row per TR/RT subject pair.
    # Entries are natural-log(Test / Reference), not raw PK measurements.
    first = np.array([[.10, -.05], [.20, .00], [-.10, .10], [.15, -.15],
                      [.00, .05], [-.05, .20], [.30, -.10], [.10, .00],
                      [.05, .15], [-.10, .10], [.20, -.05], [.00, .10]])
    estimate = analyse(first)
    decision = interim(estimate, design)
    print("Stage 1:", estimate)
    print("Decision:", decision)
    if decision.action == "continue":
        print(f"Collect {decision.additional_n} additional subjects, then use:")
        print("combined = analyse(first, second)")
        print("final(combined, design)")
    if decision.alpha is not None:
        print(f"{100*(1-2*decision.alpha):.2f}% CI:", estimate.ci(decision.alpha))


if __name__ == "__main__":
    main()
