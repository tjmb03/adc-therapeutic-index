"""Tumor PK/PD for an ADC -- efficacy driven by conjugate exposure.

A minimal tumor-growth-inhibition (TGI) model: exponential growth opposed by a
kill term proportional to conjugate concentration (the species that delivers
payload to the tumor). Repeated dosing is applied as instantaneous conjugate
increments; between doses the conjugate decays with its effective elimination
rate ``ke_eff = ke + kdec``.

The teaching point the figure makes: under-dosing shows *regrowth between doses*
even at the same cumulative dose -- so schedule, not just total dose, matters, and
you model the time course rather than a single AUC.
"""
from __future__ import annotations

import numpy as np

__all__ = ["simulate_tgi"]


def simulate_tgi(dose, t, *, ke_eff, kg, kkill, V0=200.0, dose_times=(0.0,)):
    """Tumor volume over time for a given per-administration conjugate ``dose``.

    ``kg`` tumor growth rate (1/day), ``kkill`` kill potency per unit conjugate
    per day, ``V0`` baseline volume, ``dose_times`` the dosing schedule (days).
    Integrated with a small fixed step so dosing events land exactly.
    """
    t = np.asarray(t, dtype=float)
    dt = 0.02
    n = int(round((t[-1] - t[0]) / dt))
    grid = t[0] + dt * np.arange(n + 1)
    Cc = 0.0
    V = float(V0)
    dts = sorted(dose_times)
    out_V, out_Cc = [], []
    j = 0
    for tv in grid:
        while j < len(dts) and tv >= dts[j] - dt / 2:
            Cc += dose
            j += 1
        out_V.append(V)
        out_Cc.append(Cc)
        Cc += (-ke_eff * Cc) * dt
        V += (kg * V - kkill * Cc * V) * dt
        V = max(V, 1e-6)
    out_V = np.array(out_V)
    out_Cc = np.array(out_Cc)
    # resample onto requested t
    V_t = np.interp(t, grid, out_V)
    Cc_t = np.interp(t, grid, out_Cc)
    return {"t": t, "tumor_volume": V_t, "conjugate": Cc_t}
