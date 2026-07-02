"""Multi-analyte PK behaviour and the lumped-vs-DAR identifiability result."""
import numpy as np
import pytest

from adcti import simulate_lumped, simulate_dar


T = np.linspace(0, 21, 500)
KE, KS, DAR0 = 0.10, 0.14, 4


def test_deconjugation_conserves_total_antibody():
    """With no clearance, deconjugation only moves mass between DAR bins."""
    D = simulate_dar(1000, T, ke_of_dar=lambda i: 0.0, k_single=KS,
                     DAR0=DAR0, DAR_max=8, kel_payload=6.0)
    assert np.allclose(D["total_ab"], 1000.0, rtol=1e-4)


def test_lumped_equals_dar_under_dar_independent_clearance():
    """The identifiability result: identical conjugated-payload PK."""
    L = simulate_lumped(1000, T, ke=KE, kdec=KS, DAR0=DAR0, kel_payload=6.0)
    D = simulate_dar(1000, T, ke_of_dar=lambda i: KE, k_single=KS,
                     DAR0=DAR0, DAR_max=8, kel_payload=6.0)
    rel = np.abs(L["conjugated_payload"] - D["conjugated_payload"]) / L["conjugated_payload"][0]
    assert rel.max() < 1e-4


def test_dar_matches_analytical_solution():
    """Conjugated payload = DAR0*dose*exp(-(ke+k)t); average DAR = DAR0*exp(-k t)."""
    D = simulate_dar(1000, T, ke_of_dar=lambda i: KE, k_single=KS,
                     DAR0=DAR0, DAR_max=8, kel_payload=6.0)
    cp = DAR0 * 1000 * np.exp(-(KE + KS) * T)
    assert np.allclose(D["conjugated_payload"], cp, rtol=1e-3)
    assert np.allclose(D["avg_dar"], DAR0 * np.exp(-KS * T), rtol=1e-3)


def test_dar_dependent_clearance_diverges_from_lumped():
    """The models separate only when clearance varies with DAR."""
    L = simulate_lumped(1000, T, ke=KE, kdec=KS, DAR0=DAR0, kel_payload=6.0)
    D = simulate_dar(1000, T, ke_of_dar=lambda i: KE * (1 + 0.18 * i), k_single=KS,
                     DAR0=DAR0, DAR_max=8, kel_payload=6.0)
    rel = np.abs(L["conjugated_payload"] - D["conjugated_payload"]) / L["conjugated_payload"][0]
    assert rel.max() > 0.02


def test_released_payload_rises_then_falls():
    """Free payload is formation-rate-limited: an interior Cmax."""
    L = simulate_lumped(1000, T, ke=KE, kdec=KS, DAR0=DAR0, kel_payload=6.0)
    p = L["payload"]
    i = int(np.argmax(p))
    assert 0 < i < len(p) - 1
    assert p[i] > p[0] and p[i] > p[-1]


def test_conjugate_never_exceeds_total_antibody():
    """The conjugate-vs-total-Ab gap (payload shedding) is always >= 0."""
    L = simulate_lumped(1000, T, ke=KE, kdec=KS, DAR0=DAR0, kel_payload=6.0)
    assert np.all(L["conjugate"] <= L["total_ab"] + 1e-9)
