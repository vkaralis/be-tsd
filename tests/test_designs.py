import math
import pytest
from twostage_be import Design, Estimate, interim, final


def e(n=24, gmr=1.0, cv=.2):
    return Estimate(n, math.log(gmr), math.log1p(cv**2), n-2)


def test_default_method_definitions():
    assert Design("B").assumed_gmr == .95
    assert Design("C").adjusted_alpha == .0294
    assert Design("D").assumed_gmr == .90
    assert Design("D").adjusted_alpha == .028
    assert Design("TSD-1", 150).adjusted_alpha == .028
    assert Design("TSD-2", 150).adjusted_alpha == .0294


@pytest.mark.parametrize("method", ["TSD-1", "TSD-2"])
def test_observed_gmr_gate(method):
    result = interim(e(gmr=1.3), Design(method, 150))
    assert result.action == "fail"
    assert result.reason == "gmr_outside_limits"


@pytest.mark.parametrize("method", ["B", "C", "D"])
def test_potvin_does_not_apply_observed_gmr_gate(method):
    result = interim(e(n=12, gmr=1.3, cv=.6), Design(method))
    assert result.action == "continue"
    assert result.additional_n >= 2


@pytest.mark.parametrize("method", ["TSD-1", "C", "D"])
def test_high_power_branch_stops_even_after_failure(method):
    design = Design(method, 150 if method.startswith("TSD") else None)
    passed = interim(e(n=60, cv=.1), design)
    assert passed.action == "pass"
    assert passed.alpha == .05
    failed = interim(e(n=60, gmr=1.24, cv=.1), design)
    # TSD-1 replans using 1.24; test the sufficient-power failure for fixed GMR only.
    if method != "TSD-1":
        assert failed.action == "fail"
        assert failed.reason == "sufficient_power"


@pytest.mark.parametrize("method", ["B", "TSD-2"])
def test_be_first_branch_and_power_failure(method):
    design = Design(method, 150 if method.startswith("TSD") else None)
    passed = interim(e(n=60, cv=.1), design)
    assert passed.action == "pass"
    assert passed.alpha == .0294
    if method == "B":
        failed = interim(e(n=60, gmr=1.24, cv=.1), design)
        assert failed.action == "fail"
        assert failed.reason == "sufficient_power"


def test_upper_limit_equality_is_allowed_and_excess_rejected():
    estimate = e(n=12, gmr=1.05, cv=.4)
    unlimited = interim(estimate, Design("TSD-2", 500))
    assert unlimited.action == "continue"
    total = estimate.n + unlimited.additional_n
    equal = interim(estimate, Design("TSD-2", total))
    assert equal.action == "continue"
    below = interim(estimate, Design("TSD-2", total-1))
    assert below.action == "fail"
    assert below.reason == "upper_limit"
    assert below.additional_n == 0


def test_n1_equal_or_greater_than_ul_never_recruits():
    for ul in [10, 12, 13]:
        result = interim(e(n=12, cv=.6), Design("TSD-1", ul))
        assert result.action == "fail"
        assert result.reason == "upper_limit"


def test_final_uses_adjusted_alpha_and_stage_df():
    estimate = Estimate(60, 0.0, math.log1p(.2**2), 57)
    assert final(estimate, Design("D")).action == "pass"
    assert final(estimate, Design("D")).alpha == .028
    with pytest.raises(ValueError):
        final(e(), Design("B"))
    with pytest.raises(ValueError):
        interim(estimate, Design("B"))


@pytest.mark.parametrize("kwargs", [{"method": "A"}, {"method": "TSD-1"},
    {"method": "TSD-2", "upper_limit": 150, "assumed_gmr": .95},
    {"method": "B", "assumed_gmr": 1.25}, {"method": "B", "upper_limit": -1},
    {"method": "B", "power_engine": "exact"}])
def test_invalid_design(kwargs):
    with pytest.raises(ValueError):
        Design(**kwargs)
