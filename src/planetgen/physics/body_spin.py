# planetgen/physics/body_spin.py

"""
The spin of a generated planet or moon under the torques on it (GEN.957)
========================================================================

`spin_torques` is the physics of the torques on a spinning body (the tide of the star it goes round, the tide of each moon, the thermal tide
of its air) and the flow of a spin under them; this module fills it from a body's stored fields and runs it for the age of the body's star:

* **Rheology.** A rocky body's `k2/Q` is `tidal_rheology.rock_response` at the semi-diurnal frequency of the tide it feels, so the torque depends
  on how fast the body turns (`spin_torques.rheology_dissipation`), with the mantle temperature, melt and viscosity of `interior_heat` (the
  temperature is `interior_heat.tide_temperature_k`, the same that sets the circularisation and the heat) and the solidus of its pressure. A
  gas giant keeps its fixed `k2` and `Q` (`tidal_heating`). An icy body's depends on what its water makes of it (`icy_regime`, GEN.1755, GEN.1756): the ice shell's own (`icy_tide`, GEN.960) where
  a shell sits on an ocean, solid ice's (`tidal_rheology.ice_response`) where the shell would be thicker than the body, and the rock's, with the air's thermal tide, where the surface is too warm for ice.
* **The star.** The first raiser is what the body goes round, a star for a planet and the planet for a moon, on the body's own orbit.
* **The air.** A rocky planet with an atmosphere feels the thermal tide of its star on it [R2606]: `spin_torques.ThermalTide` with the
  amplitude and relaxation time of the stored surface pressure and the starlight at the body's orbit. Where the thermal and gravitational tides
  balance away from synchronism the planet keeps a spin of its own (Venus, retrograde; the non-synchronous planets of the inner habitable zone).
* **The state.** The spin is carried from its birth spin until a synchronous state, a spin-orbit resonance whose window holds the orbit's
  eccentricity, or a balance of the tides catches it, or the age is up (`SpinResult`).

The moons of a planet act on its spin in the same sum (`torques_of` takes them) and move as they do (`moon_orbits`, GEN.959). Stand-ins: the
rows of `docs/design/research/gap-filling-plan.md` sections 6.14 (SP7 to SP10, now with a frequency-dependent torque), 6.23 and 6.24.
"""

import math

from planetgen.physics import aquifers, constants, hydrosphere, icy_tide, interior_heat, radiogenic, spin_torques, tidal_heating, tidal_rheology, tidal_spin
from planetgen.physics import eclipses, radiation_field, tide_sources


def planet_of(body):
    """The `spin_torques.Planet` of a body: its mass and radius with the moment of inertia of its kind."""
    return spin_torques.Planet(body.mass, body.radius * 1000.0, tidal_spin.moment_factor(body))


SHELL, FROZEN, ROCK = "shell", "frozen", "rock"
"""str: How an icy body (`tidal_heating.is_icy`) flexes (`icy_regime`): an ice shell over an ocean, ice all the way through, or (no ice at its surface) the rock under its water."""


def icy_regime(body, planet):
    """
    Which tide model a body's spin sees (GEN.1755, GEN.1756). `ROCK` for a body that is not icy, and for an icy one whose surface is at or above the melting point of ice (no shell: the
    water is a sea over a rocky mantle, and the mantle's rheology is what the spin tide meets, with the air's thermal tide on it, where the shell model gave 0 and the planet never despun).
    `SHELL` for a shell that sits on an ocean (`hydrosphere.ice_shell_thickness_km` under the radius, the case `icy_tide.k2_over_q` is made for). `FROZEN` for a shell thicker than the radius
    (a cold, small body: the ice is solid through, and the shell formula, which assumes an ocean under it, gave a `k2/Q` of 1 where solid ice's is a thousandth). [Stand-ins 6.23 IC6 (frozen, tier T6, tweak L1) and IC7 (warm, tier T6, tweak L0).]
    """
    if not tidal_heating.is_icy(body):
        return ROCK
    flux, top = interior_heat.icy_shell_flux_w_m2(body), interior_heat.icy_shell_top_k(body)
    thickness_m = hydrosphere.ice_shell_thickness_km(flux, top) * 1000.0
    if thickness_m <= 0.0:
        return ROCK
    return FROZEN if thickness_m >= planet.radius_m else SHELL


def frozen_mean_temperature_k(radius_m, flux_w_m2, top_k):
    """The mean of the surface temperature and the center's of a conducting body of radius `radius_m` (`T = T_top exp(F r / A)`, `k = A / T`), K: the one temperature `tidal_rheology.ice_response` is run at."""
    center_k = min(top_k * math.exp(flux_w_m2 * radius_m / hydrosphere.ICE_CONDUCTIVITY_A_W_M), aquifers.WATER_MELTING_K)
    return 0.5 * (top_k + center_k)


def dissipation_of(body, planet):
    """
    The `k2/Q` of `body` as a function of the tide's frequency, the callable `spin_torques.Torques` takes: a gas giant's or an icy body's fixed pair, and for rock the
    rheology at the temperature of `interior_heat.tide_temperature_k` with the viscosity it drew and the solidus and liquidus of the pressure in its mantle.
    """
    if getattr(body, "body_type", None) == "g":
        return spin_torques.constant_dissipation(tidal_heating.GIANT_K2, tidal_heating.GIANT_TIDAL_Q)
    regime = icy_regime(body, planet)
    if regime == SHELL:
        flux, top = interior_heat.icy_shell_flux_w_m2(body), interior_heat.icy_shell_top_k(body)
        return spin_torques.tabulate(lambda frequency: icy_tide.k2_over_q(planet.mass_kg, planet.radius_m, frequency, flux, top))
    if regime == FROZEN:
        # GEN.1755: frozen through, no ocean for the shell model to stand on: the body flexes as the ice of its own temperature does
        flux, top = interior_heat.icy_shell_flux_w_m2(body), interior_heat.icy_shell_top_k(body)
        mean_k = frozen_mean_temperature_k(planet.radius_m, flux, top)
        return spin_torques.tabulate(lambda frequency: tidal_rheology.ice_response(planet.mass_kg, planet.radius_m, frequency, mean_k).k2_over_q)
    interior = interior_heat.interior_of(body)
    if interior is None or radiogenic.content_of(body) is None:
        raise interior_heat.InteriorHeatError("a rocky body's spin needs its stored content and mantle draws")
    mantle = interior_heat.mantle_of(body.mass, body.radius, radiogenic.rock_of(body))
    return spin_torques.rheology_dissipation(planet, interior_heat.tide_temperature_k(body), interior.eta0_pa_s,
                                             interior.activation_kj_mol * 1000.0, mantle.solidus_k, mantle.liquidus_k)


def starlight_w_m2(body):
    """The light the body gets, W/m2, averaged over its orbit (the circle of the same mean flux, `orbit_elements.flux_distance_au`): 0 with no star."""
    from planetgen.physics import orbit_elements
    star = getattr(body, "star", None)
    luminosity = getattr(star, "luminosity", None)
    if not luminosity or not body.distance:
        return 0.0
    distance_au = orbit_elements.flux_distance_au(body.distance, body.host_eccentricity if body.is_moon else body.eccentricity)
    return (luminosity / (4.0 * math.pi * (distance_au * constants.AU_TO_M) ** 2) + radiation_field.companion_light_w_m2(star, distance_au)) * eclipses.kept_share(eclipses.light_factor_of(body))     # GEN.1441, GEN.1484


def thermal_tide_of(body, star_raiser):
    """The `spin_torques.ThermalTide` of a rocky planet's air under its star, or `None` for a body with no air, a moon (it is held by its planet) or a body that
    is not rock."""
    if getattr(body, "body_type", None) != "t" or getattr(body, "is_moon", False) or icy_regime(body, planet_of(body)) != ROCK:      # an icy body with a shell or ice through it keeps none; one too warm for ice has the air's tide (GEN.1756)
        return None
    pressure_pa = getattr(body, "atmospheric_pressure", None)
    flux = starlight_w_m2(body)
    if not pressure_pa or pressure_pa <= 0.0 or flux <= 0.0:
        return None
    return spin_torques.ThermalTide(star_raiser, spin_torques.thermal_pressure_amplitude_pa(pressure_pa, flux),
                                    spin_torques.thermal_relaxation_rate(pressure_pa, flux))


def torques_of(body, host_mass_kg, semi_major_axis_m, moons=()):
    """The `spin_torques.Torques` on `body`: the tide of what it orbits (first), the tide of each moon in `moons` (`spin_torques.Raiser`s), the tide of every other body of the system that pulls on it (`tide_sources.spin_raisers`, GEN.1502) and the thermal tide of its air."""
    planet = planet_of(body)
    star = spin_torques.Raiser(host_mass_kg, semi_major_axis_m, spin_torques.orbital_rate(host_mass_kg, semi_major_axis_m))
    return spin_torques.Torques(planet, (star, *moons, *tide_sources.spin_raisers(body)), dissipation_of(body, planet), thermal_tide_of(body, star))


def evolve(body, host_mass_kg, semi_major_axis_m, spin_rate, age_s, eccentricity=0.0, moons=(), capture_units=spin_torques.ALWAYS_CAPTURED):
    """
    A body born spinning at `spin_rate` (rad/s) after `age_s` of the torques on it: the `spin_torques.SpinResult`, with the rate in rad/s (negative for a
    retrograde spin), the ratio of the spin to the orbital rate if a synchronous state or a resonance has caught it, and the kind of state.
    """
    return spin_torques.evolve(torques_of(body, host_mass_kg, semi_major_axis_m, moons), spin_rate, age_s, eccentricity, capture_units)
