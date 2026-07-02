"""Multi-analyte ADC PK and the lumped-vs-DAR identifiability result.

    python examples/multianalyte_pk.py
"""
import numpy as np

from adcti import simulate_lumped, simulate_dar

T = np.linspace(0, 21, 400)
KE, KS, DAR0 = 0.10, 0.14, 4


def main():
    L = simulate_lumped(1000, T, ke=KE, kdec=KS, DAR0=DAR0, kel_payload=6.0)
    Dind = simulate_dar(1000, T, ke_of_dar=lambda i: KE, k_single=KS,
                        DAR0=DAR0, DAR_max=8, kel_payload=6.0)
    Ddep = simulate_dar(1000, T, ke_of_dar=lambda i: KE * (1 + 0.20 * i), k_single=KS,
                        DAR0=DAR0, DAR_max=8, kel_payload=6.0)

    rel_ind = np.abs(L["conjugated_payload"] - Dind["conjugated_payload"]).max() / L["conjugated_payload"][0]
    rel_dep = np.abs(L["conjugated_payload"] - Ddep["conjugated_payload"]).max() / L["conjugated_payload"][0]

    print("Multi-analyte ADC PK (lumped model), t = 0 -> 21 d")
    print(f"  conjugate    Cmax {L['conjugate'][0]:.0f} -> {L['conjugate'][-1]:.1f} nM")
    print(f"  total Ab          {L['total_ab'][0]:.0f} -> {L['total_ab'][-1]:.1f} nM")
    print(f"  free payload peak {L['payload'].max():.2f} nM (formation-rate-limited)")
    print(f"  average DAR  {Dind['avg_dar'][0]:.2f} -> {Dind['avg_dar'][-1]:.3f}")
    print()
    print("Lumped vs DAR-resolved (total conjugated payload):")
    print(f"  DAR-INDEPENDENT clearance : max rel diff {rel_ind:.2e}   -> identical / unidentifiable")
    print(f"  DAR-DEPENDENT   clearance : max rel diff {rel_dep:.3f}     -> divergent / identifiable")


if __name__ == "__main__":
    main()
