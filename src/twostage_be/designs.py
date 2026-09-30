"""Explicit decision rules for the five supported two-stage methods."""
from dataclasses import dataclass
import math

from .statistics import Estimate, planning_power, required_n

METHODS = ("TSD-1", "TSD-2", "B", "C", "D")


@dataclass(frozen=True)
class Design:
    method: str
    upper_limit: int | None = None
    assumed_gmr: float | None = None
    target_power: float = 0.8
    power_engine: str | None = None

    def __post_init__(self):
        if self.method not in METHODS:
            raise ValueError(f"method must be one of {METHODS}")
        if self.upper_limit is not None and (
            type(self.upper_limit) is not int or self.upper_limit < 1
        ):
            raise ValueError("upper_limit must be a positive integer or None")
        if self.method.startswith("TSD"):
            if self.upper_limit is None:
                raise ValueError("TSD-1/2 require a predefined upper_limit")
            if self.assumed_gmr is not None:
                raise ValueError("TSD-1/2 use observed GMR; do not supply assumed_gmr")
        elif self.assumed_gmr is None:
            object.__setattr__(self, "assumed_gmr", 0.90 if self.method == "D" else 0.95)
        if self.assumed_gmr is not None and (
            not math.isfinite(self.assumed_gmr) or not 0.8 < self.assumed_gmr < 1.25
        ):
            raise ValueError("assumed_gmr must lie strictly between 0.8 and 1.25")
        if not math.isfinite(self.target_power) or not 0 < self.target_power < 1:
            raise ValueError("target_power must be between 0 and 1")
        if self.power_engine is None:
            object.__setattr__(self, "power_engine",
                               "legacy" if self.method.startswith("TSD") else "potvin")
        if self.power_engine not in {"legacy", "potvin"}:
            raise ValueError("power_engine must be legacy or potvin")

    @property
    def adjusted_alpha(self):
        return 0.028 if self.method in {"TSD-1", "D"} else 0.0294


@dataclass(frozen=True)
class Decision:
    action: str  # pass, fail, continue
    reason: str
    alpha: float | None
    additional_n: int = 0
    estimated_power: float | None = None


def interim(estimate: Estimate, design: Design) -> Decision:
    """Decide whether stage 1 passes, fails, or recruits a second stage."""
    if estimate.df != estimate.n - 2:
        raise ValueError("Interim estimate must represent a single balanced stage")
    adaptive_gmr = design.method.startswith("TSD")
    if adaptive_gmr and not math.log(0.8) <= estimate.log_gmr <= math.log(1.25):
        return Decision("fail", "gmr_outside_limits", None)
    gmr = estimate.gmr if adaptive_gmr else design.assumed_gmr
    alpha = design.adjusted_alpha
    engine = design.power_engine
    if design.method in {"TSD-1", "C", "D"}:
        power = planning_power(estimate.n, estimate.mse, gmr, 0.05, engine)
        if power >= design.target_power:
            passed = estimate.equivalent(0.05)
            return Decision("pass" if passed else "fail", "sufficient_power", 0.05,
                            estimated_power=power)
        if estimate.equivalent(alpha):
            return Decision("pass", "stage1_be", alpha, estimated_power=power)
    else:
        if estimate.equivalent(alpha):
            return Decision("pass", "stage1_be", alpha)
        power = planning_power(estimate.n, estimate.mse, gmr, alpha, engine)
        if power >= design.target_power:
            return Decision("fail", "sufficient_power", alpha, estimated_power=power)
    if design.upper_limit is not None and estimate.n + 2 > design.upper_limit:
        return Decision("fail", "upper_limit", alpha, estimated_power=power)
    needed = required_n(estimate.mse, gmr, alpha, design.target_power,
                        engine, design.upper_limit)
    if needed is None:
        return Decision("fail", "upper_limit", alpha, estimated_power=power)
    total = max(needed, estimate.n + 2)
    if design.upper_limit is not None and total > design.upper_limit:
        return Decision("fail", "upper_limit", alpha, estimated_power=power)
    return Decision("continue", "stage2_required", alpha, total - estimate.n, power)


def final(estimate: Estimate, design: Design) -> Decision:
    """Assess the combined, stage-adjusted analysis once at the adjusted alpha."""
    if estimate.df != estimate.n - 3:
        raise ValueError("Final estimate must represent two balanced stages")
    alpha = design.adjusted_alpha
    passed = estimate.equivalent(alpha)
    return Decision("pass" if passed else "fail", "combined_be" if passed else "combined_ci",
                    alpha)
