"""Multi-analyte pharmacokinetics for an antibody-drug conjugate (ADC).

An ADC is never one PK curve. At minimum three species move together:

* **conjugate** -- antibody still carrying >=1 payload; drives tumor delivery / efficacy;
* **total antibody** -- conjugated + fully deconjugated antibody; the long-lived backbone;
* **released payload** -- free cytotoxin, cleared fast on its own volume; the safety driver.

Deconjugation sheds payload over time, so the average DAR falls and the conjugate
curve converges toward total antibody; higher-DAR species may also clear faster.

Two model granularities are provided:

* :func:`simulate_lumped` -- a single conjugate pool with an average DAR (3 analytes);
* :func:`simulate_dar` -- a DAR0..DAR_max cascade with per-DAR clearance.

**Identifiability result (pinned in the tests).** Under *DAR-independent* clearance,
the total conjugated-payload PK of the DAR-resolved model is identical to the lumped
model's ``conjugate x DAR`` -- both decay as ``exp(-(ke + k_deconj) t)`` and the average
DAR falls as ``DAR0 * exp(-k_deconj t)``. The two models **diverge only when clearance
varies with DAR.** So the DAR-resolved model is only identifiable with DAR-distribution
data (or demonstrably DAR-dependent clearance); otherwise the lumped model suffices and
the extra parameters are unsupported.
"""
from __future__ import annotations

import numpy as np
from scipy.integrate import odeint

__all__ = ["simulate_lumped", "simulate_dar"]


def simulate_lumped(dose, t, *, ke, kdec, DAR0, kel_payload,
                    k12=0.0, k21=0.0, payload_vr=10.0):
    """Lumped three-analyte ADC PK (conjugate / total antibody / released payload).

    Parameters are first-order rate constants (1/day); ``dose`` is the initial
    conjugate concentration; ``payload_vr`` is the payload's volume ratio
    (larger => lower free-payload concentration). Returns a dict of arrays.
    """
    def rhs(y, _t):
        Cc, Cp, Dc, Dp, Pl = y
        dCc = -(ke + kdec + k12) * Cc + k21 * Cp
        dCp = k12 * Cc - k21 * Cp
        dDc = kdec * Cc - (ke + k12) * Dc + k21 * Dp
        dDp = k12 * Dc - k21 * Dp
        dPl = kdec * Cc * DAR0 / payload_vr - kel_payload * Pl
        return [dCc, dCp, dDc, dDp, dPl]

    Y = odeint(rhs, [dose, 0.0, 0.0, 0.0, 0.0], t)
    Cc, Cp, Dc, Dp, Pl = Y.T
    return {
        "t": t,
        "conjugate": Cc,
        "total_ab": Cc + Dc,
        "payload": Pl,
        "conjugated_payload": Cc * DAR0,   # payload still attached (single-pool avg DAR)
    }


def simulate_dar(dose, t, *, ke_of_dar, k_single, DAR0, DAR_max, kel_payload,
                 payload_vr=10.0):
    """DAR-resolved ADC PK: a DAR0..DAR_max deconjugation cascade.

    Antibody starts monodisperse at ``DAR0``. Each attached payload deconjugates
    independently at ``k_single`` (so a DAR-i species leaves for DAR-(i-1) at rate
    ``i * k_single``). ``ke_of_dar`` is a callable ``i -> clearance`` (constant for
    DAR-independent clearance, increasing for DAR-dependent clearance).
    """
    N = DAR_max
    idx = np.arange(N + 1)
    ke_i = np.array([float(ke_of_dar(i)) for i in idx])

    def rhs(y, _t):
        dar = y[:N + 1]
        Pl = y[N + 1]
        d = np.zeros(N + 2)
        deconj = idx * k_single * dar          # DAR i -> i-1 flux
        d[:N + 1] -= ke_i * dar + deconj
        d[:N] += deconj[1:]                     # DAR i-1 gains from DAR i
        d[N + 1] = deconj.sum() / payload_vr - kel_payload * Pl
        return d

    y0 = np.zeros(N + 2)
    y0[DAR0] = dose
    Y = odeint(rhs, y0, t)
    dar = Y[:, :N + 1]
    Pl = Y[:, N + 1]
    total_ab = dar.sum(1)
    conj_payload = (dar * idx).sum(1)
    conjugate = dar[:, 1:].sum(1)
    avg_dar = np.divide(conj_payload, total_ab,
                        out=np.zeros_like(conj_payload), where=total_ab > 0)
    return {
        "t": t,
        "dar": dar,
        "total_ab": total_ab,
        "conjugate": conjugate,
        "conjugated_payload": conj_payload,
        "avg_dar": avg_dar,
        "payload": Pl,
    }
