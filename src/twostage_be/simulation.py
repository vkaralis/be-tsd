"""Reproducible Monte Carlo studies for one bioequivalence endpoint."""
from collections import Counter
from dataclasses import asdict, dataclass
import math

import numpy as np
from scipy.stats import norm

from .designs import Design, interim, final
from .statistics import analyse


@dataclass(frozen=True)
class Trial:
    passed: bool
    stage: int
    total_n: int
    reason: str
    alpha: float | None
    observed_gmr: float
    observed_cv: float


def _population(gmr, cv):
    if not math.isfinite(gmr) or gmr <= 0 or not math.isfinite(cv) or cv < 0:
        raise ValueError("gmr must be positive; cv must be nonnegative; both finite")


def generate_stage(n: int, gmr: float, cv: float, rng: np.random.Generator,
                   period_effect: float = 0.0) -> np.ndarray:
    """Generate paired log differences in TR/RT columns.

    The residual variance of a difference is 2*log(1+CVw**2).
    period_effect is the period-1 minus period-2 difference on the log scale.
    """
    _population(gmr, cv)
    if type(n) is not int or n < 2 or n % 2:
        raise ValueError("n must be a positive even integer >=2")
    if not math.isfinite(period_effect):
        raise ValueError("period_effect must be finite")
    return rng.normal(math.log(gmr), math.sqrt(2 * math.log1p(cv**2)), (n // 2, 2)) + np.array(
        [period_effect, -period_effect])


def simulate_trial(design: Design, n1: int, gmr: float, cv: float,
                   rng: np.random.Generator) -> Trial:
    if type(n1) is not int or n1 < 4 or n1 % 2:
        raise ValueError("n1 must be an even integer >=4")
    first = generate_stage(n1, gmr, cv, rng)
    estimate = analyse(first)
    decision = interim(estimate, design)
    stage = 1
    if decision.action == "continue":
        second = generate_stage(decision.additional_n, gmr, cv, rng)
        estimate = analyse(first, second)
        decision = final(estimate, design)
        stage = 2
    return Trial(decision.action == "pass", stage, estimate.n, decision.reason,
                 decision.alpha, estimate.gmr, estimate.cv)


def wilson_interval(successes: int, total: int) -> tuple[float, float]:
    """95% Wilson interval for the Monte Carlo acceptance probability."""
    if total <= 0 or not 0 <= successes <= total:
        raise ValueError("Invalid binomial counts")
    z = float(norm.ppf(0.975))
    p = successes / total
    den = 1 + z*z / total
    center = (p + z*z / (2 * total)) / den
    half = z / den * math.sqrt(p*(1-p)/total + z*z/(4*total*total))
    return max(0.0, center - half), min(1.0, center + half)


def simulate(design: Design, n1: int = 24, gmr: float = 0.95, cv: float = 0.2,
             trials: int = 1000, seed: int = 2026) -> dict:
    """Run a scenario; proportions are fractions, never percentages.

    mean_total_n includes all trials, including those rejected before recruitment.
    At true GMR=0.8 or 1.25, acceptance estimates boundary type I error.
    """
    if type(trials) is not int or trials < 1:
        raise ValueError("trials must be a positive integer")
    if type(n1) is not int or n1 < 4 or n1 % 2:
        raise ValueError("n1 must be an even integer >=4")
    _population(gmr, cv)
    rng = np.random.default_rng(seed)
    successes = stage2 = sample_sum = stage1_successes = stage2_successes = 0
    reasons = Counter()
    for _ in range(trials):
        trial = simulate_trial(design, n1, gmr, cv, rng)
        successes += trial.passed
        stage2 += trial.stage == 2
        stage1_successes += trial.passed and trial.stage == 1
        stage2_successes += trial.passed and trial.stage == 2
        sample_sum += trial.total_n
        reasons[trial.reason] += 1
    probability = successes / trials
    low, high = wilson_interval(successes, trials)
    return {
        **asdict(design), "n1": n1, "true_gmr": gmr, "true_cv": cv,
        "trials": trials, "seed": seed, "accepted": successes,
        "acceptance_probability": probability,
        "mc_standard_error": math.sqrt(probability*(1-probability)/trials),
        "mc_ci_low": low, "mc_ci_high": high,
        "stage2_probability": stage2 / trials,
        "stage1_acceptance_probability": stage1_successes / trials,
        "stage2_acceptance_probability": stage2_successes / trials,
        "mean_total_n": sample_sum / trials,
        "reasons": dict(sorted(reasons.items())),
    }
