"""Screen the ADC panel for therapeutic window and OBD, ranked by index.

    python examples/dose_selection.py

Also regenerates figures/adc_panel.png.
"""
import os

import pandas as pd

from adcti import screen_adc_panel
from adcti.figures import fig_adc_panel

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CSV = os.path.join(ROOT, "data", "adc_candidates.csv")


def main():
    df = screen_adc_panel(CSV)
    show = df[[
        "name", "start_dose", "MED", "MTD", "OBD",
        "efficacy_at_MTD", "therapeutic_index", "dosable", "verdict",
    ]].copy()
    for c in ["start_dose", "MED", "MTD", "OBD", "efficacy_at_MTD", "therapeutic_index"]:
        show[c] = show[c].round(3)

    pd.set_option("display.width", 240)
    pd.set_option("display.max_columns", 30)
    print("\nADC dose-selection screen — ranked by therapeutic index\n")
    print(show.to_string(index=False))

    out = os.path.join(ROOT, "figures", "adc_panel.png")
    fig_adc_panel(out)
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
