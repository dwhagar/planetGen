# planetgen/physics/heavy_gas_loss.py

"""
The loss of the heavy gases of an air through the star's XUV history (GEN.1558, a prebuild: nothing calls it yet)
=================================================================================================================

The stepped pressure caps of the retention regime (`atmosphere_recipes` row `d6`, read by `surface_air._retention_cap_pa`) stand in for a loss that
`thermosphere.kept_shares` does not follow: its Jeans flux is nil for CO2, N2 and Ar, and its ion channel (`ion_escape`) is off in the air stage
because the field (`magnetopause_rp`) is drawn after the air. This module is the loss term that replaces the caps: a pure function of the body, the
star's XUV history and the field, with no draw and no state.

1. **What the XUV could lift.** `lifted_mass_kg` is the sum over the history of `thermosphere.energy_limited_mass_loss_kg_s` (the energy-limited rate
   at an efficiency of 1, at 1.05 planet radii as `kept_shares` does for its ion channel), on the same log-spaced steps and the same XUV decay
   (`thermosphere.xuv_ratio`, `_mixed_history` for the companions). It is an upper bound: every absorbed photon lifting gas.
2. **What of it is the air's.** Heavy gases leave only by being dragged, or as ions the wind sweeps up, so a share `efficiency` of the bound goes
   to the air, over the surface the field leaves open (`ion_escape.open_share`). The air is a column that nothing resupplies here, so the share
   kept is `exp(-efficiency * open * lifted / air_mass)`.

The measured spread of `lifted / air` is 10 to 3000 for a Mars, a Venus, an Earth and a Moon-sized body at their present pressures (a bound, not
a loss: it is 850 for the Earth, which kept its air). A single efficiency separates the Earth (open share 0.05 at ten radii) from Mars (open all
over) only in a narrow window, near 0.005 to 0.015. That window is the stand-in the build has to register (GEN.725) and calibrate; this module takes
the efficiency as an argument and invents none.
"""

import math

from planetgen import tuning
from planetgen.physics import constants, ion_escape, thermosphere

HISTORY_STEPS = thermosphere.HISTORY_STEPS
"""int: Log-spaced steps of the history (the same as `thermosphere.kept_shares`)."""
ABSORBING_RADIUS_OVER_PLANET = 1.05
"""float: The radius that absorbs the XUV, in planet radii (the one `kept_shares` uses for the ion channel)."""


def lifted_mass_kg(mass_kg, radius_m, xuv_flux_now, age_gy, host_mass_sol=None, luminosity_w=None, log_lx_lbol=None, companions=()):
    """
    The most mass the XUV of the body's history could lift out of its well by `age_gy`, kg (the energy-limited bound; see the module).

    Args:
        mass_kg, radius_m (float): The body's.
        xuv_flux_now (float): `thermosphere.xuv_flux_earth` of the body, in the Earth's present.
        age_gy (float): The body's age.
        host_mass_sol, luminosity_w (float, optional): The star's, for the shape of its decay; without them the flux is held at its present value.
        log_lx_lbol, companions: As in `thermosphere.kept_shares`.
    """
    if age_gy <= 0.0 or xuv_flux_now <= 0.0:
        return 0.0
    start_gy = min(tuning.RUNAWAY_HISTORY_START_YR / 1e9, 0.5 * age_gy)
    ratio_step = (age_gy / start_gy) ** (1.0 / HISTORY_STEPS)
    t0, total = start_gy, 0.0
    for _ in range(HISTORY_STEPS):
        t1 = t0 * ratio_step
        mid = math.sqrt(t0 * t1)
        history = thermosphere.xuv_ratio(host_mass_sol, luminosity_w, age_gy, mid, log_lx_lbol) if host_mass_sol and luminosity_w else 1.0
        if companions:
            history = thermosphere._mixed_history(history, xuv_flux_now, companions, age_gy, mid)
        flux = min(xuv_flux_now * history, thermosphere.XUV_FLUX_MAX_EARTH) * constants.EARTH_XUV_FLUX_W_M2
        rate = thermosphere.energy_limited_mass_loss_kg_s(flux, ABSORBING_RADIUS_OVER_PLANET * radius_m, mass_kg, efficiency=1.0)
        total += rate * (t1 - t0) * thermosphere._GY_S
        t0 = t1
    return total


def kept_share(air_mass_kg, lifted_kg, efficiency, magnetopause_rp=None):
    """
    The share of an air of `air_mass_kg` that stays, 0..1: `exp(-efficiency * open_share * lifted_kg / air_mass_kg)`, `open_share` that of
    `ion_escape` (all of the surface for no field, `None` or 0). A body with no air keeps nothing.
    """
    if air_mass_kg <= 0.0:
        return 0.0
    if efficiency < 0.0:
        raise ValueError("the efficiency of the heavy-gas loss cannot be negative")
    return math.exp(-min(efficiency * ion_escape.open_share(magnetopause_rp) * lifted_kg / air_mass_kg, 700.0))
