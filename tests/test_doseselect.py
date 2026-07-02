"""The dose-selection window and the toxicology-paradigm decision rule."""
import numpy as np
import pytest

from adcti import efficacy, toxicity, select_dose


def test_med_and_mtd_closed_form():
    r = select_dose("t", ED50=1.0, f_payload=0.05, TD50=0.5, hnstd_hed=6.0)
    # MED at 50% efficacy = ED50
    assert r.MED == pytest.approx(1.0, rel=1e-6)
    # MTD: free-payload exposure at 30% DLT = TD50*(0.3/0.7)**0.5, then / f_payload
    dp = 0.5 * (0.3 / 0.7) ** 0.5
    assert r.MTD == pytest.approx(dp / 0.05, rel=1e-6)


def test_therapeutic_index_is_mtd_over_med():
    r = select_dose("t", ED50=1.0, f_payload=0.05, TD50=0.5, hnstd_hed=6.0)
    assert r.therapeutic_index == pytest.approx(r.MTD / r.MED, rel=1e-9)


def test_start_dose_is_hed_over_safety_factor():
    r = select_dose("t", ED50=1.0, f_payload=0.05, TD50=0.5, hnstd_hed=6.0, safety_factor=6.0)
    assert r.start_dose == pytest.approx(1.0)


def test_unstable_linker_lowers_mtd():
    """More free payload per dose (higher f_payload) tightens the ceiling."""
    stable = select_dose("s", ED50=1.0, f_payload=0.05, TD50=0.5, hnstd_hed=6.0)
    unstable = select_dose("u", ED50=1.0, f_payload=0.15, TD50=0.5, hnstd_hed=6.0)
    assert unstable.MTD < stable.MTD
    assert unstable.therapeutic_index < stable.therapeutic_index


def test_no_window_is_no_go():
    """MED above MTD => no dose is both efficacious and safe."""
    r = select_dose("t", ED50=5.0, f_payload=0.15, TD50=0.3, hnstd_hed=6.0)
    assert r.dosable is False
    assert r.MED > r.MTD
    assert r.therapeutic_index < 1.0


def test_obd_at_or_below_mtd():
    for ED50 in [0.3, 1.0, 3.0]:
        r = select_dose("t", ED50=ED50, f_payload=0.05, TD50=0.5, hnstd_hed=6.0)
        assert r.OBD <= r.MTD + 1e-9


def test_efficacy_and_toxicity_monotone():
    d = np.logspace(-2, 2, 400)
    assert np.all(np.diff(efficacy(d, ED50=1.0)) > 0)
    assert np.all(np.diff(toxicity(d, f_payload=0.05, TD50=0.5)) > 0)
    assert efficacy(1e3, 1.0) > 0.99 and efficacy(1e-3, 1.0) < 0.01
