"""Regenerate every figure from the model. Run ``python -m adcti.figures``.

Nothing here is hand-drawn: each panel is computed from :mod:`adcti.pk`,
:mod:`adcti.pd` and :mod:`adcti.doseselect`. The two ADC schematics
(``adc_lumped_scheme.png``, ``adc_dar_scheme.png``) are committed static assets.
"""
from __future__ import annotations

import os

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from .pk import simulate_lumped, simulate_dar
from .pd import simulate_tgi
from .doseselect import efficacy, toxicity, select_dose, screen_adc_panel

TEAL, CORAL, AMBER, PURPLE, GREY = "#0E8C6B", "#D63B30", "#C8860B", "#6A4FB0", "#8A8A82"
BAND = "#E9E6DD"

plt.rcParams.update({
    "font.size": 11, "axes.edgecolor": "#8A8A82", "axes.linewidth": 1.0,
    "axes.spines.top": False, "axes.spines.right": False,
    "xtick.color": "#5F5E5A", "ytick.color": "#5F5E5A",
    "axes.labelcolor": "#2C2C2A", "figure.dpi": 150,
})

HERE = os.path.dirname(os.path.abspath(__file__))
FIGDIR = os.path.normpath(os.path.join(HERE, "..", "..", "figures"))
DATA = os.path.normpath(os.path.join(HERE, "..", "..", "data", "adc_candidates.csv"))

# reference ADC for the dose-selection map (distinct start/MED/OBD/MTD)
REF = dict(ED50=1.0, f_payload=0.04, TD50=0.5, hnstd_hed=3.0)


def fig_multianalyte_pk(path):
    t = np.linspace(0, 28, 600)
    L = simulate_lumped(1000, t, ke=0.10, kdec=0.14, DAR0=4, kel_payload=6.0,
                        k12=0.25, k21=0.18)
    fig, ax = plt.subplots(figsize=(7.4, 3.8))
    ax.semilogy(t, L["total_ab"], color=TEAL, lw=2.6, label="total antibody")
    ax.semilogy(t, L["conjugate"], color=PURPLE, lw=2.6, label="conjugate (ADC)")
    ax.semilogy(t, np.clip(L["payload"], 1e-3, None), color=CORAL, lw=2.4, label="released payload")
    ig = np.argmin(np.abs(t - 18))
    ax.annotate("payload shed\n(DAR \u2193)", xy=(18, L["conjugate"][ig]),
                xytext=(19, L["total_ab"][ig] * 0.9), fontsize=9, color="#5F5E5A",
                bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.85))
    ax.set_xlabel("time after dose (days)")
    ax.set_ylabel("concentration (nM, log)")
    ax.set_ylim(0.1, 2000)
    ax.legend(loc="upper right", frameon=False, fontsize=9.5, ncol=3, bbox_to_anchor=(1.0, 1.15))
    fig.tight_layout(); fig.savefig(path, bbox_inches="tight"); plt.close(fig)


def fig_lumped_vs_dar(path):
    t = np.linspace(0, 21, 400)
    ke, ks, DAR0 = 0.10, 0.14, 4
    L = simulate_lumped(1000, t, ke=ke, kdec=ks, DAR0=DAR0, kel_payload=6.0)
    Dind = simulate_dar(1000, t, ke_of_dar=lambda i: ke, k_single=ks,
                        DAR0=DAR0, DAR_max=8, kel_payload=6.0)
    Ddep = simulate_dar(1000, t, ke_of_dar=lambda i: ke * (1 + 0.20 * i), k_single=ks,
                        DAR0=DAR0, DAR_max=8, kel_payload=6.0)
    fig, ax = plt.subplots(figsize=(7.4, 3.9))
    ax.semilogy(t, L["conjugated_payload"], color=TEAL, lw=2.8, label="lumped model")
    ax.semilogy(t[::14], Dind["conjugated_payload"][::14], "o", color=TEAL, ms=6,
                mfc="white", mew=1.6, label="DAR-resolved, DAR-independent CL")
    ax.semilogy(t, Ddep["conjugated_payload"], color=CORAL, lw=2.4, ls="--",
                label="DAR-resolved, DAR-dependent CL")
    ax.annotate("identical\n(unidentifiable)", xy=(6, L["conjugated_payload"][114]),
                xytext=(2.0, 250), fontsize=9.5, color="#0C5E47",
                bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.85),
                arrowprops=dict(arrowstyle="-", color=GREY, lw=0.7, ls=":"))
    ax.annotate("diverges only when\nCL varies with DAR", xy=(15, Ddep["conjugated_payload"][285]),
                xytext=(11.5, 60), fontsize=9.5, color="#9A3520",
                bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.85),
                arrowprops=dict(arrowstyle="-", color=GREY, lw=0.7, ls=":"))
    ax.set_xlabel("time after dose (days)")
    ax.set_ylabel("total conjugated payload (nM, log)")
    ax.legend(loc="upper right", frameon=True, facecolor="white", framealpha=0.9,
              edgecolor="none", fontsize=9)
    fig.tight_layout(); fig.savefig(path, bbox_inches="tight"); plt.close(fig)


def fig_tgi(path):
    t = np.linspace(0, 35, 700)
    doses = [(0, GREY, "--", "vehicle"), (100, AMBER, "-", "1 mg/kg"),
             (300, PURPLE, "-", "3 mg/kg"), (1000, TEAL, "-", "10 mg/kg")]
    fig, ax = plt.subplots(figsize=(7.4, 3.8))
    for d, col, ls, lab in doses:
        r = simulate_tgi(d, t, ke_eff=0.24, kg=0.07, kkill=0.00075, V0=200.0,
                         dose_times=(0, 7, 14, 21))
        ax.plot(t, r["tumor_volume"], color=col, ls=ls, lw=2.6, label=lab)
    for dt in (0, 7, 14, 21):
        ax.plot([dt, dt], [0, 60], color="#B23A2E", lw=2.2)
    ax.text(1.2, 2350, "red ticks = doses (weekly \u00D74)", color="#9A3520", fontsize=9)
    ax.set_xlabel("time (days)")
    ax.set_ylabel("tumor volume (mm\u00B3)")
    ax.set_ylim(0, 2500)
    ax.legend(loc="upper left", frameon=False, fontsize=9.5, ncol=4, bbox_to_anchor=(0.0, 1.14))
    fig.tight_layout(); fig.savefig(path, bbox_inches="tight"); plt.close(fig)


def fig_dose_selection_map(path):
    r = select_dose("ref", **REF)
    d = np.logspace(-2, 2, 600)
    E = efficacy(d, REF["ED50"]) * 100
    Tx = toxicity(d, REF["f_payload"], REF["TD50"]) * 100
    fig, ax = plt.subplots(figsize=(7.6, 4.0))
    ax.axvspan(r.MED, r.MTD, color=BAND, alpha=0.9, lw=0)
    ax.semilogx(d, E, color=TEAL, lw=3, label="efficacy \u2014 conjugate exposure")
    ax.semilogx(d, Tx, color=CORAL, lw=3, label="DLT rate \u2014 free payload")
    ax.axhline(50, color=TEAL, ls="--", lw=0.8, alpha=0.5)
    ax.axhline(30, color=CORAL, ls="--", lw=0.8, alpha=0.5)
    ax.text(1.05e-2, 52, "efficacy target", color="#0C5E47", fontsize=9)
    ax.text(1.05e-2, 31.5, "max acceptable DLT", color="#9A3520", fontsize=9)
    marks = [(r.start_dose, AMBER, "start (HED)", 4), (r.MED, TEAL, "MED", 13),
             (r.OBD, PURPLE, "OBD", 4), (r.MTD, CORAL, "MTD", 13)]
    for x, c, lab, yy in marks:
        ax.axvline(x, color=c, ls="--", lw=1.1)
        ax.text(x, yy, lab, color=c, fontsize=8.5, ha="center",
                bbox=dict(boxstyle="round,pad=0.1", fc="white", ec="none", alpha=0.85))
    ax.text((r.MED * r.MTD) ** 0.5, 66, "therapeutic\nwindow", ha="center", fontsize=9.5,
            color="#5a5a55",
            bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.85))
    ax.set_xlabel("dose (mg/kg, log)")
    ax.set_ylabel("response / probability (%)")
    ax.set_ylim(0, 100)
    ax.legend(loc="center left", frameon=True, facecolor="white", framealpha=0.9,
              edgecolor="none", fontsize=9.5, bbox_to_anchor=(0.02, 0.62))
    fig.text(0.5, -0.02, "Start & ceiling from payload toxicology; floor from conjugate efficacy; "
             "OBD chosen in-window \u2014 no MABEL.", ha="center", fontsize=9, color="#5F5E5A")
    fig.tight_layout(); fig.savefig(path, bbox_inches="tight"); plt.close(fig)


def fig_adc_panel(path):
    df = screen_adc_panel(DATA).iloc[::-1].reset_index(drop=True)  # widest TI on top
    fig, ax = plt.subplots(figsize=(8.6, 4.6))
    for i, r in df.iterrows():
        lo, hi = min(r.MED, r.MTD), max(r.MED, r.MTD)
        if not r.dosable:
            col = CORAL
        elif r.therapeutic_index >= 5:
            col = TEAL
        else:
            col = AMBER
        ax.plot([lo, hi], [i, i], color=col, lw=6, solid_capstyle="round",
                alpha=0.5 if not r.dosable else 1.0)
        ax.plot([r.MED], [i], "o", color=TEAL if r.dosable else CORAL, ms=8, zorder=5,
                label="MED (efficacy floor)" if i == 0 else None)
        ax.plot([r.MTD, r.MTD], [i - 0.26, i + 0.26], color=CORAL, lw=2.2,
                label="MTD (payload ceiling)" if i == 0 else None)
        if r.dosable:
            ax.plot([r.OBD], [i], "D", color=PURPLE, ms=6, zorder=6,
                    label="OBD / RP2D" if i == 0 else None)
        ax.plot([r.start_dose], [i], "|", color=AMBER, ms=14, mew=2.4,
                label="start dose (HED)" if i == 0 else None)
        tag = "GO" if (r.dosable and r.therapeutic_index >= 5) else ("GO*" if r.dosable else "NO-GO")
        ax.text(1.02, i, f"{tag}  (TI {r.therapeutic_index:.1f})", va="center", ha="left",
                fontsize=9, color=col, transform=ax.get_yaxis_transform())

    ax.set_xscale("log")
    ax.set_yticks(range(len(df)))
    ax.set_yticklabels(df["name"])
    ax.set_ylim(-0.6, len(df) - 0.4)
    ax.set_xlim(0.1, 30)
    ax.set_xlabel("dose (mg/kg, log)")
    ax.set_title("ADC dose-selection screen \u2014 therapeutic window (MED \u2192 MTD), ranked by index",
                 fontsize=11.5, color="#2C2C2A", loc="left", pad=10)
    ax.legend(loc="upper center", bbox_to_anchor=(0.42, -0.13), ncol=4,
              frameon=False, fontsize=8.5)
    fig.subplots_adjust(right=0.76, bottom=0.2)
    fig.savefig(path, bbox_inches="tight"); plt.close(fig)


def main():
    os.makedirs(FIGDIR, exist_ok=True)
    fig_multianalyte_pk(os.path.join(FIGDIR, "multianalyte_pk.png"))
    fig_lumped_vs_dar(os.path.join(FIGDIR, "lumped_vs_dar.png"))
    fig_tgi(os.path.join(FIGDIR, "tgi.png"))
    fig_dose_selection_map(os.path.join(FIGDIR, "dose_selection_map.png"))
    fig_adc_panel(os.path.join(FIGDIR, "adc_panel.png"))
    print("wrote figures to", FIGDIR)


if __name__ == "__main__":
    main()
