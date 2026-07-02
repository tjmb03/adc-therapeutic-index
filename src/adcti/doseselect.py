"""Exposure-driven dose selection for a cytotoxic ADC -- the toxicology paradigm.

For a cytotoxic ADC the dose-limiting toxicity is **off-target free payload**, so
dose selection follows the toxicology / exposure-response paradigm, **not MABEL**:

* the **start dose** is toxicology-anchored -- an animal HNSTD (or STD10) carried to a
  human-equivalent dose and divided by a safety factor (default 6, the ICH S9 HNSTD
  approach);
* the therapeutic window is bounded **below** by the **conjugate** exposure-response
  (minimum efficacious dose, MED) and **above** by the **free-payload** exposure-response
  (maximum tolerated dose, MTD);
* the **OBD / RP2D** is chosen *inside* that window by exposure-response (Project Optimus),
  not at the MTD.

Conjugate target occupancy lives on the **efficacy** side -- it sets the active dose, not
the start dose. (The one place conjugate occupancy re-enters as a *safety* input is
on-target/off-tumor binding, e.g. CEACAM5 on normal GI epithelium -- handled outside this
model.) Contrast with an immune-agonist bispecific, where engagement itself is the hazard
and MABEL governs; see the companion repo ``bispecific-fih-dosability``.

The efficacy curve is driven by conjugate exposure (approximately linear in dose here);
the toxicity curve by free-payload exposure ``dose * f_payload``, where ``f_payload`` is the
fraction of dose appearing as free payload (linker instability / deconjugation).
"""
from __future__ import annotations

from dataclasses import dataclass, asdict

import numpy as np
import pandas as pd

__all__ = ["DoseWindow", "efficacy", "toxicity", "select_dose", "screen_adc_panel"]


def efficacy(dose, ED50, h=1.5):
    """Fractional efficacy vs dose (conjugate-driven Emax)."""
    dose = np.asarray(dose, dtype=float)
    e = dose ** h / (dose ** h + ED50 ** h)
    return float(e) if e.ndim == 0 else e


def toxicity(dose, f_payload, TD50, h=2.0):
    """Fractional DLT probability vs dose (free-payload-driven Emax)."""
    dose = np.asarray(dose, dtype=float)
    dp = dose * f_payload
    tox = dp ** h / (dp ** h + TD50 ** h)
    return float(tox) if tox.ndim == 0 else tox


def _dose_for_effect(level, EC50, h):
    return EC50 * (level / (1.0 - level)) ** (1.0 / h)


@dataclass
class DoseWindow:
    name: str
    ED50: float
    f_payload: float
    TD50: float
    hnstd_hed: float
    start_dose: float          # HED / safety factor (toxicology-anchored)
    MED: float                 # min efficacious dose (conjugate exposure-response)
    MTD: float                 # max tolerated dose (free-payload exposure-response)
    OBD: float                 # optimal biological dose / RP2D, chosen in-window
    efficacy_at_MTD: float
    therapeutic_index: float   # MTD / MED
    dosable: bool
    verdict: str

    def as_row(self) -> dict:
        return asdict(self)


def _verdict(dosable, TI, eff_at_mtd):
    if not dosable:
        return "NO-GO: no therapeutic window (MED above MTD)"
    if eff_at_mtd < 0.5:
        return "GO*: window is narrow / tox-limited below strong efficacy"
    if TI >= 5:
        return "GO: wide therapeutic index"
    return "GO: workable window"


def select_dose(name, ED50, f_payload, TD50, hnstd_hed,
                eff_target=0.5, dlt_max=0.30, h_eff=1.5, h_tox=2.0,
                safety_factor=6.0, obd_efficacy=0.90):
    """Select the FIH window and OBD for one ADC; see module docstring."""
    MED = _dose_for_effect(eff_target, ED50, h_eff)

    dp_at_dlt = _dose_for_effect(dlt_max, TD50, h_tox)   # free-payload exposure at DLT
    MTD = dp_at_dlt / f_payload                          # convert to dose

    start = hnstd_hed / safety_factor

    # OBD: efficacy-plateau dose, capped at the MTD (Project Optimus: not the MTD by default)
    dose_at_plateau = _dose_for_effect(obd_efficacy, ED50, h_eff)
    OBD = min(dose_at_plateau, MTD)

    eff_at_mtd = float(efficacy(MTD, ED50, h_eff))
    TI = MTD / MED
    dosable = (MED < MTD) and (start < MTD)

    return DoseWindow(
        name=name, ED50=ED50, f_payload=f_payload, TD50=TD50, hnstd_hed=hnstd_hed,
        start_dose=start, MED=MED, MTD=MTD, OBD=OBD, efficacy_at_MTD=eff_at_mtd,
        therapeutic_index=TI, dosable=dosable, verdict=_verdict(dosable, TI, eff_at_mtd),
    )


def screen_adc_panel(candidates, **kw) -> pd.DataFrame:
    """Screen a panel of ADCs, ranked by therapeutic index (widest first).

    ``candidates`` is a CSV path (columns
    ``name,ED50_mgkg,f_payload,TD50_payload,HNSTD_HED_mgkg[,note]``) or a DataFrame.
    """
    if isinstance(candidates, (str, bytes)) or hasattr(candidates, "__fspath__"):
        df = pd.read_csv(candidates)
    else:
        df = pd.DataFrame(candidates)

    rows = [
        select_dose(
            name=r["name"], ED50=r["ED50_mgkg"], f_payload=r["f_payload"],
            TD50=r["TD50_payload"], hnstd_hed=r["HNSTD_HED_mgkg"], **kw,
        )
        for _, r in df.iterrows()
    ]
    out = pd.DataFrame([r.as_row() for r in rows])
    if "note" in df.columns:
        out = out.merge(df[["name", "note"]], on="name", how="left")
    return out.sort_values("therapeutic_index", ascending=False).reset_index(drop=True)
