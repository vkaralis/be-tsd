import numpy as np
import pytest
from twostage_be import Design, simulate, generate_stage, analyse
from twostage_be.simulation import wilson_interval


@pytest.mark.parametrize("method", ["TSD-1", "TSD-2", "B", "C", "D"])
def test_reproducible_bounded_simulation(method):
    design = Design(method, 62)
    a = simulate(design, n1=12, cv=.4, trials=60, seed=718)
    b = simulate(design, n1=12, cv=.4, trials=60, seed=718)
    assert a == b
    assert 12 <= a["mean_total_n"] <= 62
    assert sum(a["reasons"].values()) == 60
    assert a["acceptance_probability"] == pytest.approx(
        a["stage1_acceptance_probability"] + a["stage2_acceptance_probability"])
    assert a["mc_ci_low"] <= a["acceptance_probability"] <= a["mc_ci_high"]


def test_population_and_period_effect():
    data = generate_stage(20000, 1.05, .3, np.random.default_rng(820), period_effect=.4)
    estimate = analyse(data)
    assert estimate.log_gmr == pytest.approx(np.log(1.05), abs=.01)
    assert estimate.mse == pytest.approx(np.log1p(.3**2), rel=.04)


def test_binomial_interval_extremes():
    assert wilson_interval(0, 100)[0] == pytest.approx(0, abs=1e-15)
    assert wilson_interval(100, 100)[1] == pytest.approx(1)


@pytest.mark.parametrize("kwargs", [{"trials": 0}, {"n1": 13}, {"cv": -1},
                                    {"gmr": float("nan")}])
def test_invalid_scenario(kwargs):
    with pytest.raises(ValueError):
        simulate(Design("B"), **kwargs)
