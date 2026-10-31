# planetgen/physics/body_spin.py

"""
The spin of a generated planet or moon under the torques on it (GEN.957)
========================================================================

`spin_torques` is the physics of the torques on a spinning body (the tide of the star it goes round, the tide of each moon, the thermal tide
of its air) and the flow of a spin under them; this module fills it from a body's stored fields and runs it for the age of the body's star:

* **Rheology.** A rocky body's `k2/Q` is `tidal_rheology.rock_response` at the semi-diurnal frequency of the tide it feels, so the torque depends
  on how fast the body turns (`spin_torques.rheology_dissipation`), with the mantle temperature, melt and viscosity of `interior_heat` (the
  temperature is `interior_heat.tide_temperature_k`, the same that sets the circularisation and the heat) and the solidus of its pressure. A
  gas giant and an icy body keep their fixed `k2` and `Q` (`tidal_heating`, GEN.960 gives ice its own).
* **The star.** The first raiser is what the body goes round, a star for a planet and the planet for a moon, on the body's own orbit.
* **The air.** A rocky planet with an atmosphere feels the thermal tide of its star on it (Leconte et al. 2015): `spin_torques.ThermalTide` with the
  amplitude and relaxation time of the stored surface pressure and the starlight at the body's orbit. Where the thermal and gravitational tides
  balance away from synchronism the planet keeps a spin of its own (Venus, retrograde; the non-synchronous planets of the inner habitable zone).
* **The state.** The spin is carried from its birth spin until a synchronous state, a spin-orbit resonance whose window holds the orbit's
  eccentricity, or a balance of the tides catches it, or the age is up (`SpinResult`).

The moons of a planet act on its spin in the same sum (`torques_of` takes them) and move as they do (`moon_orbits`, GEN.959). Stand-ins: the
rows of `docs/design/research/gap-filling-plan.md` sections 6.14 (SP7 to SP10, now with a frequency-dependent torque), 6.23 and 6.24.
"""

import math

from planetgen.physics import constants, icy_tide, interior_heat, radiogenic, spin_torques, tidal_heating, tidal_spin


def planet_of(body):
    """The `spin_torques.Planet` of a body: its mass and radius with the moment of inertia of its kind."""
    return spin_torques.Planet(body.mass, body.radius * 1000.0, tidal_spin.moment_factor(body))


def dissipation_of(body, planet):
    """
    The `k2/Q` of `body` as a function of the tide's frequency, the callable `spin_torques.Torques` takes: a gas giant's or an icy body's fixed pair, and for rock the
    rheology at the temperature of `interior_heat.tide_temperature_k` with the viscosity it drew and the solidus and liquidus of the pressure in its mantle.
    """
    if getattr(body, "body_type", None) == "g":
        return spin_torques.constant_dissipation(tidal_heating.GIANT_K2, tidal_heating.GIANT_TIDAL_Q)
    if tidal_heating.is_icy(body):
        flux, top = interior_heat.icy_shell_flux_w_m2(body), interior_heat.icy_shell_top_k(body)
        return spin_torques.tabulate(lambda frequency: icy_tide.k2_over_q(planet.mass_kg, planet.radius_m, frequency, flux, top))
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
    return luminosity / (4.0 * math.pi * (distance_au * constants.AU_TO_M) ** 2)


def thermal_tide_of(body, star_raiser):
    """The `spin_torques.ThermalTide` of a rocky planet's air under its star, or `None` for a body with no air, a moon (it is held by its planet) or a body that
    is not rock."""
    if getattr(body, "body_type", None) != "t" or getattr(body, "is_moon", False) or tidal_heating.is_icy(body):
        return None
    pressure_pa = getattr(body, "atmospheric_pressure", None)
    flux = starlight_w_m2(body)
    if not pressure_pa or pressure_pa <= 0.0 or flux <= 0.0:
        return None
    return spin_torques.ThermalTide(star_raiser, spin_torques.thermal_pressure_amplitude_pa(pressure_pa, flux),
                                    spin_torques.thermal_relaxation_rate(pressure_pa, flux))


def torques_of(body, host_mass_kg, semi_major_axis_m, moons=()):
    """The `spin_torques.Torques` on `body`: the tide of what it orbits (first), the tide of each moon in `moons` (`spin_torques.Raiser`s) and the thermal tide of its air."""
    planet = planet_of(body)
    star = spin_torques.Raiser(host_mass_kg, semi_major_axis_m, spin_torques.orbital_rate(host_mass_kg, semi_major_axis_m))
    return spin_torques.Torques(planet, (star, *moons), dissipation_of(body, planet), thermal_tide_of(body, star))


def evolve(body, host_mass_kg, semi_major_axis_m, spin_rate, age_s, eccentricity=0.0, moons=(), capture_units=spin_torques.ALWAYS_CAPTURED):
    """
    A body born spinning at `spin_rate` (rad/s) after `age_s` of the torques on it: the `spin_torques.SpinResult`, with the rate in rad/s (negative for a
    retrograde spin), the ratio of the spin to the orbital rate if a synchronous state or a resonance has caught it, and the kind of state.
    """
    return spin_torques.evolve(torques_of(body, host_mass_kg, semi_major_axis_m, moons), spin_rate, age_s, eccentricity, capture_units)
