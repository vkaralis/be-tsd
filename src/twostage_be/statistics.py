"""Statistics for complete, balanced two-sequence crossover stages.

Differences are ln(Test) - ln(Reference), one value per subject.
"""
from dataclasses import dataclass
from functools import lru_cache
import math

import numpy as np
from scipy.stats import norm, t


@lru_cache(maxsize=8192)
def critical_t(alpha: float, df: int) -> float:
    return float(t.ppf(1 - alpha, df))


@dataclass(frozen=True)
class Estimate:
    n: int
    log_gmr: float
    mse: float
    df: int

    def __post_init__(self):
        if self.n < 4 or self.n % 2 or not 0 < self.df <= self.n - 2:
            raise ValueError("Invalid sample size or residual degrees of freedom")
        if not math.isfinite(self.log_gmr) or not math.isfinite(self.mse) or self.mse < 0:
            raise ValueError("Estimates must be finite and MSE nonnegative")

    @property
    def gmr(self) -> float:
        return math.exp(self.log_gmr)

    @property
    def cv(self) -> float:
        return math.sqrt(math.expm1(self.mse))

    def log_ci(self, alpha: float) -> tuple[float, float]:
        _alpha(alpha)
        half = critical_t(alpha, self.df) * math.sqrt(2 * self.mse / self.n)
        return self.log_gmr - half, self.log_gmr + half

    def ci(self, alpha: float) -> tuple[float, float]:
        return tuple(math.exp(x) for x in self.log_ci(alpha))

    def equivalent(self, alpha: float) -> bool:
        lo, hi = self.log_ci(alpha)
        return lo >= math.log(0.8) and hi <= math.log(1.25)


def _alpha(alpha):
    if not math.isfinite(alpha) or not 0 < alpha < 0.5:
        raise ValueError("alpha must be between 0 and 0.5")


def analyse(*stages) -> Estimate:
    """Fit treatment + period(stage) to paired differences.

Each stage is an (n/2, 2) array; columns are TR and RT sequences.
This is the paired-difference reduction of the fixed subject ANOVA.
The common treatment coefficient is the overall mean in this balanced case.
Stage-specific sequence contrasts represent period effects. Stage and
subject intercepts cancel on differencing. Two-stage residual df = N - 3.
"""
    if not 1 <= len(stages) <= 2:
        raise ValueError("Supply one or two stages")
    arrays = [np.asarray(stage, dtype=float) for stage in stages]
    for a in arrays:
        if a.ndim != 2 or a.shape[1] != 2 or len(a) < 1 or not np.isfinite(a).all():
            raise ValueError("Each stage must be a finite (n/2, 2) array")
    n = sum(a.size for a in arrays)
    df = n - len(arrays) - 1
    if df <= 0:
        raise ValueError("Insufficient residual degrees of freedom")
    mean = sum(float(a.sum()) for a in arrays) / n
    ss = 0.0
    for a in arrays:
        # Remove the stage-specific period contrast, retaining the common mean.
        contrast = (float(a[:, 0].mean()) - float(a[:, 1].mean())) / 2
        residual = a - mean - np.array([contrast, -contrast])
        ss += float(np.sum(residual**2)) / 2
    return Estimate(n, mean, ss / df, df)


def planning_power(n: int, mse: float, gmr: float, alpha: float,
                   engine: str = "potvin", df: int | None = None) -> float:
    """Published planning approximations, not exact unconditional TOST power.

    potvin: central-t CDF difference, Potvin (2008), section 2.3.
    legacy: nearest-boundary Julious correction used in supplied MATLAB.
    Planning uses n-2 df by default; actual tests use the fitted residual df.
    """
    _alpha(alpha)
    if n < 4 or n % 2 or not math.isfinite(mse) or mse < 0:
        raise ValueError("n must be even and >=4; mse must be finite and >=0")
    if not math.isfinite(gmr) or gmr <= 0:
        raise ValueError("gmr must be finite and positive")
    if engine not in {"potvin", "legacy"}:
        raise ValueError("Unknown power engine")
    df = n - 2 if df is None else df
    if df <= 0:
        raise ValueError("df must be positive")
    delta = math.log(gmr)
    distance = min(delta - math.log(0.8), math.log(1.25) - delta)
    if distance <= 0:
        return 0.0
    if mse == 0:
        return 1.0
    q = critical_t(alpha, df)
    if engine == "legacy":
        z = -q + math.sqrt(max(0, n - q / 2) * distance**2 / (2 * mse))
        return float(t.cdf(z, df))
    se = math.sqrt(2 * mse / n)
    hi = (math.log(1.25) - delta) / se - q
    lo = (math.log(0.8) - delta) / se + q
    return float(max(0, t.cdf(hi, df) - t.cdf(lo, df)))


def required_n(mse: float, gmr: float, alpha: float, target: float = 0.8,
               engine: str = "potvin", max_n: int | None = None) -> int | None:
    """Required total size; None means the requirement exceeds max_n.

    Potvin uses a monotone search for the smallest even n attaining target.
    Legacy follows the supplied Julious iterative formula and MATLAB rounding.
    The upper limit rejects a requirement; it never truncates recruitment.
    """
    planning_power(4, mse, gmr, alpha, engine)
    if not math.isfinite(target) or not 0 < target < 1:
        raise ValueError("target must be between 0 and 1")
    if max_n is not None and (not isinstance(max_n, int) or max_n < 1):
        raise ValueError("max_n must be a positive integer")
    distance = min(math.log(gmr / 0.8), math.log(1.25 / gmr))
    if distance <= 0:
        return None
    if engine == "legacy":
        beta_quantile = target if gmr != 1 else (1 + target) / 2
        za, zb = norm.ppf(1 - alpha), norm.ppf(beta_quantile)
        x = math.floor(2 * mse * (za + zb)**2 / distance**2 + za / 2) + 1
        if x <= 2:
            final = x + 2
        else:
            for _ in range(100):
                df = math.floor(x + 0.5) - 2
                qa, qb = critical_t(alpha, df), float(t.ppf(beta_quantile, df))
                final = 2 * mse * (qa + qb)**2 / distance**2 + qa / 2
                # One-sided tolerance deliberately follows the historical code.
                if final - x < 0.5:
                    break
                x = final
            else:
                raise RuntimeError("Legacy sample-size iteration did not converge")
        n = max(4, 2 * math.ceil(math.floor(final + 0.5) / 2))
        return None if max_n is not None and n > max_n else n
    if planning_power(4, mse, gmr, alpha, engine) >= target:
        return None if max_n is not None and max_n < 4 else 4
    lo, hi = 2, 4  # half sample sizes
    cap = max_n // 2 if max_n is not None else None
    if cap is not None and cap < 2:
        return None
    while planning_power(2 * hi, mse, gmr, alpha, engine) < target:
        lo = hi
        if cap is not None and hi >= cap:
            return None
        hi = min(2 * hi, cap) if cap is not None else 2 * hi
    if cap is not None and hi > cap:
        hi = cap
        if planning_power(2 * hi, mse, gmr, alpha, engine) < target:
            return None
    while hi - lo > 1:
        mid = (hi + lo) // 2
        if planning_power(2 * mid, mse, gmr, alpha, engine) >= target:
            hi = mid
        else:
            lo = mid
    return 2 * hi
