import math
import numpy as np
import pytest
from scipy.stats import t

from twostage_be import analyse, planning_power, required_n


@pytest.mark.parametrize("sizes", [(12,), (12, 2), (24, 18)])
def test_difference_analysis_matches_full_fixed_subject_anova(sizes):
    """Independent full observation-level GLM, including one-pair stage 2."""
    rng = np.random.default_rng(781)
    responses, treatments, periods, subjects, differences = [], [], [], [], []
    offset = 0
    for stage, n in enumerate(sizes):
        y = rng.normal(size=(n, 2))
        seq = np.tile([1, -1], n // 2)
        # First column Test, second Reference; sequence determines period.
        responses.extend(y.reshape(-1))
        treatments.extend(np.tile([1, 0], n))
        for j, sign in enumerate(seq):
            subjects.extend([offset + j, offset + j])
            periods.extend([(stage, int(sign == 1)), (stage, int(sign != 1))])
        differences.append((y[:, 0] - y[:, 1]).reshape(n // 2, 2))
        offset += n
    # Full subject intercepts absorb sequence/stage and their nested subjects.
    x = np.column_stack([
        np.eye(offset)[subjects], np.asarray(treatments),
        *[np.array([int(s == stage and p == 1) for s, p in periods])
          for stage in range(len(sizes))],
    ])
    y = np.asarray(responses)
    coefficients, _, rank, _ = np.linalg.lstsq(x, y, rcond=None)
    residual = y - x @ coefficients
    df = len(y) - rank
    mse = float(residual @ residual) / df
    estimate = analyse(*differences)
    assert estimate.df == df
    assert estimate.log_gmr == pytest.approx(coefficients[offset])
    assert estimate.mse == pytest.approx(mse)
    treatment_variance = mse * np.linalg.inv(x.T @ x)[offset, offset]
    assert treatment_variance == pytest.approx(2 * estimate.mse / estimate.n)


def test_hand_computed_stage():
    e = analyse([[0.1, -0.1], [0.3, 0.1]])
    assert e.log_gmr == pytest.approx(0.1)
    assert e.mse == pytest.approx(0.01)
    assert e.df == 2
    assert e.cv == pytest.approx(math.sqrt(math.expm1(0.01)))
    half = t.ppf(0.95, 2) * math.sqrt(0.005)
    assert e.log_ci(0.05) == pytest.approx((0.1 - half, 0.1 + half))


@pytest.mark.parametrize("gmr", [0.9, 0.95, 1.0, 1.1])
def test_potvin_size_minimal_and_ratio_reciprocal(gmr):
    mse = math.log1p(0.3**2)
    n = required_n(mse, gmr, 0.0294)
    assert n % 2 == 0
    assert planning_power(n, mse, gmr, 0.0294) >= 0.8
    if n > 4:
        assert planning_power(n - 2, mse, gmr, 0.0294) < 0.8
    assert required_n(mse, 1/gmr, 0.0294) == n
    assert required_n(mse, gmr, 0.0294, max_n=n) == n
    assert required_n(mse, gmr, 0.0294, max_n=n - 1) is None


def test_central_t_formula_and_legacy_correction():
    n, mse, gmr, alpha = 24, math.log1p(0.2**2), 0.95, 0.0294
    q = t.ppf(1-alpha, n-2)
    se = math.sqrt(2*mse/n)
    expected = t.cdf(math.log(1.25/gmr)/se-q, n-2) - t.cdf(
        math.log(0.8/gmr)/se+q, n-2)
    assert planning_power(n, mse, gmr, alpha) == pytest.approx(expected)
    expected_legacy = t.cdf(-q + math.sqrt((n-q/2)*math.log(gmr/0.8)**2/(2*mse)), n-2)
    assert planning_power(n, mse, gmr, alpha, "legacy") == pytest.approx(expected_legacy)


def test_legacy_size_numeric_fixture():
    # Julious iteration in the supplied N2.m: CVw=.2, GMR=.95, alpha=.028.
    assert required_n(math.log1p(.2**2), .95, .028, engine="legacy") == 24


@pytest.mark.parametrize("data", [[], [[1, 2, 3]], [[float("nan"), 1]], [[1, 2]]])
def test_invalid_data(data):
    with pytest.raises(ValueError):
        analyse(data)


def test_boundary_and_zero_variance():
    assert planning_power(24, .1, 1.25, .0294) == 0
    assert required_n(.1, 1.25, .0294) is None
    assert planning_power(24, 0, .95, .0294) == 1
    assert required_n(0, .95, .0294) == 4
