# planetgen/physics/mixture_escape.py

"""
Thermal escape of the light gases from a small body's mixture (GEN.693, rebuilt on the thermosphere model by GEN.738)
====================================================================================================================

A small warm world cannot keep hydrogen or helium: the gas outgassed at the surface leaves from the exobase faster than the
volcanoes supply it. The filter works on the mole fractions `mixture_draw.draw_mixture` takes from the families that outgas
(`rocky`, `venus`, `mars`, `abiotic-oxygen`); the giants, the envelopes, the lava vapour (supplied by a
magma ocean faster than it leaves), the ice-supported air and the Titan-like air are not filtered.

The share of each species that stays is the one `thermosphere.kept_shares` finds: the exobase temperature of the star's XUV and the
carbon dioxide of the air, the Jeans flux from that exobase (never more than the air can carry up by diffusion, or the XUV can lift out of
the well), summed over the star's XUV history against the species' column. The model, its sources and its limits are in
docs/design/composition/atmospheric-escape-model.md.

A species with less than `LOST_BELOW` of it left is gone; what stays is renormalised; if less than `LOST_TOTAL` of the gas stays the body
is airless. Both are numerical cutoffs, not physics: the first keeps a trace that renormalising would promote to the main gas out of the
mixture, the second keeps what is left of a mostly escaped air from being called an air. Pressure is not changed: the total comes from
the pressure stage, which already carries the loss of the atmosphere (the cosmic shoreline of rows d).
The gas that leaves is taken out of the moles of the draw (`MixtureDraw.moles`), so the pressure the supply law gives (`surface_air`) is
the pressure of what stays.
"""

import dataclasses

from planetgen.physics import constants, thermosphere
from planetgen.physics.mixture_base import MIN_FRACTION, MixtureError

FAMILIES = ("rocky", "venus", "mars", "abiotic-oxygen")
"""tuple: The families that are filtered. Not `titan` (its hydrogen trace is the steady state of photochemical production and escape,
measured by Niemann et al. 2010, which a loss-only filter would remove) and not `frost` (an air held by the vapour pressure of its ice,
which the ice resupplies faster than it leaves)."""

LOST_BELOW = 1e-3
"""float: A species of which less than this share stays is gone (a numerical cutoff: not left as a trace that renormalising would promote to the main gas)."""

LOST_TOTAL = 1e-2
"""float: A body whose gas has less than this share left after the loss is airless (a numerical cutoff: what stays of a mostly escaped air is not an air)."""


def body_inputs(subject):
    """
    `(mass_kg, radius_m, pressure_pa, age_gy)` of the subject. The radius is `radius_km`, else the one the escape speed and the mass give
    (`R = 2 G M / v^2`); a body with neither is refused by name.
    """
    age = subject.get("age_gy")
    v_esc = subject.get("v_esc_kms")
    mass_earth = subject.get("mass_earth")
    if age is None or (v_esc is None and subject.get("radius_km") is None):
        raise MixtureError("the escape filter needs the quantities v_esc_kms and age_gy")
    if mass_earth is None:
        raise MixtureError("the escape filter needs the quantity mass_earth")
    mass = mass_earth * constants.EARTH_MASS_TO_KG
    radius = subject.get("radius_km")
    radius_m = radius * 1e3 if radius else 2.0 * constants.G * mass / (v_esc * 1e3) ** 2
    pressure_bar = subject.get("pressure_bar")
    pressure_pa = (pressure_bar if pressure_bar else 1.0) * 1e5
    return mass, radius_m, pressure_pa, age


def apply(subject, result):
    """
    `result` (a `MixtureDraw`) without the gases the body loses.

    Args:
        subject (Subject): With `age_gy`, `mass_earth`, `teq_k` (else `ts_k`) and `v_esc_kms` or `radius_km`; `pressure_bar` (else 1 bar) is the
            column the loss is compared to; `xuv_flux_earth` or the host's mass, luminosity and the insolation give the XUV (see
            `thermosphere.xuv_flux_earth`; a body with no star gets none, and its thermosphere stays at the mesopause temperature).
        result (MixtureDraw): The family's mixture.

    Returns:
        MixtureDraw: The mixture renormalised, with the share each species kept in `constants["escape_kept"]`, the exobase temperature
        now in `constants["exobase_k"]` and the `moles` that stay (the moles of each species times the share it kept); an empty mixture
        (airless) if nothing stays.
    """
    if result.family not in FAMILIES or not result.mole_fractions:
        return result
    mass, radius, pressure, age = body_inputs(subject)
    teq = subject.get("teq_k") or subject.get("ts_k")
    if not teq:
        raise MixtureError("the escape filter needs the equilibrium temperature (teq_k) or the surface temperature")
    flux = thermosphere.xuv_flux_earth(subject)
    mass_sol, luminosity = subject.get("host_mass_solar"), subject.get("luminosity_solar")
    kept, exobase = thermosphere.kept_shares(result.mole_fractions, mass, radius, teq, pressure, age, flux, mass_sol,
                                             luminosity * constants.SOLAR_LUMINOSITY if luminosity else None,
                                             magnetopause_rp=subject.get("magnetopause_rp"), log_lx_lbol=subject.get("log_lx_lbol"))
    left = {formula: fraction * kept[formula] for formula, fraction in result.mole_fractions.items() if kept[formula] >= LOST_BELOW}
    left = {formula: value for formula, value in left.items() if value > MIN_FRACTION}
    total = sum(left.values())
    if total < LOST_TOTAL:
        total = 0.0
    moles = None if result.moles is None else result.moles * total
    constants_out = {**result.constants, "exobase_k": exobase, "escape_kept": {f: round(k, 6) for f, k in kept.items() if k < 0.999999}}
    if total <= 0.0:
        return dataclasses.replace(result, mole_fractions={}, family="airless", constants=constants_out)
    return dataclasses.replace(result, mole_fractions={formula: value / total for formula, value in left.items()}, constants=constants_out,
                               moles=moles)
