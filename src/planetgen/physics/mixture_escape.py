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

from planetgen.physics import constants, thermosphere, unit_registry as units
from planetgen.physics.mixture_base import MIN_FRACTION, MixtureError

_mearth_to_kg = units.converter("Mearth", "kg")

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
    mass = _mearth_to_kg(mass_earth)
    radius = subject.get("radius_km")
    radius_m = radius * 1e3 if radius else 2.0 * constants.G * mass / (v_esc * 1e3) ** 2
    pressure_bar = subject.get("pressure_bar")
    pressure_pa = (pressure_bar if pressure_bar else 1.0) * 1e5
    return mass, radius_m, pressure_pa, age


def kept_in_column(subject, fractions, pressure_pa=None):
    """
    `(kept, exobase_k)`: the share each species of an air of these mole `fractions` keeps (`thermosphere.kept_shares`) on `subject`'s body. The column the
    loss is compared to is `pressure_pa`, else the subject's `pressure_bar` (1 bar without one).
    """
    mass, radius, pressure, age = body_inputs(subject)
    if pressure_pa:
        pressure = pressure_pa
    teq = subject.get("teq_k") or subject.get("ts_k")
    if not teq:
        raise MixtureError("the escape filter needs the equilibrium temperature (teq_k) or the surface temperature")
    flux = thermosphere.xuv_flux_earth(subject)
    mass_sol, luminosity = subject.get("host_mass_solar"), subject.get("luminosity_solar")
    return thermosphere.kept_shares(fractions, mass, radius, teq, pressure, age, flux, mass_sol,
                                    luminosity * constants.SOLAR_LUMINOSITY if luminosity else None,
                                    magnetopause_rp=subject.get("magnetopause_rp"), log_lx_lbol=subject.get("log_lx_lbol"))


def survivors(fractions, kept):
    """
    `(left, total)`: the mole `fractions` times the share each species `kept`, without the species of which less than `LOST_BELOW` stays or that falls under
    `MIN_FRACTION`, and `total` the share of the air that stays (0 when it is under `LOST_TOTAL`: no air).
    """
    left = {formula: fraction * kept[formula] for formula, fraction in fractions.items() if kept[formula] >= LOST_BELOW}
    left = {formula: value for formula, value in left.items() if value > MIN_FRACTION}
    total = sum(left.values())
    return left, (0.0 if total < LOST_TOTAL else total)


def apply(subject, result, column_bar=None):
    """
    `result` (a `MixtureDraw`) without the gases the body loses.

    Args:
        subject (Subject): With `age_gy`, `mass_earth`, `teq_k` (else `ts_k`) and `v_esc_kms` or `radius_km`; `pressure_bar` (else 1 bar) is the
            column the loss is compared to; `xuv_flux_earth` or the host's mass, luminosity and the insolation give the XUV (see
            `thermosphere.xuv_flux_earth`; a body with no star gets none, and its thermosphere stays at the mesopause temperature).
        result (MixtureDraw): The family's mixture.
        column_bar (float, optional): The pressure the loss is compared to when it is not the subject's (`kept_in_column`).

    Returns:
        MixtureDraw: The mixture renormalised, with the share each species kept in `constants["escape_kept"]`, the exobase temperature
        now in `constants["exobase_k"]` and the `moles` that stay (the moles of each species times the share it kept); an empty mixture
        (airless) if nothing stays.
    """
    if result.family not in FAMILIES or not result.mole_fractions:
        return result
    kept, exobase = kept_in_column(subject, result.mole_fractions, column_bar * 1e5 if column_bar else None)
    left, total = survivors(result.mole_fractions, kept)
    moles = None if result.moles is None else result.moles * total
    constants_out = {**result.constants, "exobase_k": exobase, "escape_kept": {f: round(k, 6) for f, k in kept.items() if k < 0.999999}}
    if total <= 0.0:
        return dataclasses.replace(result, mole_fractions={}, family="airless", constants=constants_out)
    return dataclasses.replace(result, mole_fractions={formula: value / total for formula, value in left.items()}, constants=constants_out,
                               moles=moles)
