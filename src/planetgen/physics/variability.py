# planetgen/physics/variability.py

"""
The kind of variable a star is, from what the star stores (GEN.944, the cheap flags of its variable-star part).

The star sheet shows a "Variability" row for the stars that vary in a way their phase explains. Nothing is stored: the kind follows from the
star's `state` (GEN.786), temperature and mass, which the sheet reads when it is made.

| Star | Variability |
|---|---|
| state `cepheid` (a giant in the instability strip, `stellar_evolution.classical_cepheid`) | `classical Cepheid` |
| state `asymptotic-giant` cooler than `MIRA_MAX_TEMPERATURE_K` | `Mira or semiregular` |
| class `PMS` under `HERBIG_MIN_MASS_SOL` | `T Tauri` |
| class `PMS` from `HERBIG_MIN_MASS_SOL` to `HERBIG_MAX_MASS_SOL` | `Herbig Ae/Be` |

Not here: carbon and S stars (whether an asymptotic-giant star has run enough thermal pulses to turn carbon-rich is not modelled, and the
composition code only says it *can*), Ap and Am stars (they need a rotation and a magnetic field the generator gives a star only at the
end), and the long-period variables below the asymptotic-giant branch (the semiregulars of the red-giant branch).
"""

MIRA_MAX_TEMPERATURE_K = 3500.0
"""float: An asymptotic-giant star cooler than this pulsates as a Mira or semiregular variable. [Stand-in VR1, tier T6, tweak L1: Miras run from about
2,500 to 3,500 K; the hotter asymptotic-giant stars are on the way to their thermal-pulse phase and do not vary this way.]"""
HERBIG_MIN_MASS_SOL = 2.0
"""float: A pre-main-sequence star of this mass or more is a Herbig Ae/Be star, under it a T Tauri star. [Stand-in VR2, tier T6, tweak L0: the
usual class boundary, the intermediate-mass analogues of the T Tauri stars.]"""
HERBIG_MAX_MASS_SOL = 8.0
"""float: A pre-main-sequence star of this mass or more has no optical pre-main-sequence phase (it reaches the main sequence still accreting), so
it is neither. [Stand-in VR3, tier T6, tweak L0: the upper limit of the Herbig Ae/Be class.]"""

CEPHEID = "classical Cepheid"
MIRA = "Mira or semiregular"
T_TAURI = "T Tauri"
HERBIG = "Herbig Ae/Be"
KINDS = (CEPHEID, MIRA, T_TAURI, HERBIG)


def kind_of(state, yerkes_class, temperature_k, mass_sol):
    """The kind of variable a star of this `state` (`stellar_state.STATES`), Yerkes class, effective temperature (K) and mass (solar masses) is,
    one of `KINDS`, or None for a star that does not vary in a way its phase gives."""
    if state == "cepheid":
        return CEPHEID
    if state == "asymptotic-giant" and temperature_k is not None and temperature_k < MIRA_MAX_TEMPERATURE_K:
        return MIRA
    if yerkes_class == "PMS" and mass_sol is not None:
        if mass_sol < HERBIG_MIN_MASS_SOL:
            return T_TAURI
        if mass_sol < HERBIG_MAX_MASS_SOL:
            return HERBIG
    return None
