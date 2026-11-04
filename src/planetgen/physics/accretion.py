# planetgen/physics/accretion.py

"""
The X-ray luminosity of an accreting neutron star or black hole (GEN.272, docs/design/systems/multi-star-systems.md 13.4,
step 16), derived from a pair's masses, radii and orbit and never stored: `L = eta Mdot c^2`, capped at the Eddington
limit of the accretor.

Three ways the gas arrives, named by `mode`:

* `wind`: a hot massive star's wind, captured by the Bondi-Hoyle rate `(G M)^2 Mdot_wind / (a^2 v_wind v_rel^3)` with
  `v_rel^2 = v_wind^2 + v_orbit^2` (the supergiant high-mass X-ray binary). The wind loses mass at a rate that scales with
  the star's luminosity (`tuning.WIND_MASS_LOSS_MSUN_YR_AT_1E5_LSUN`) and blows at `tuning.WIND_SPEED_OVER_ESCAPE` times
  its escape speed.
* `disc`: the slow decretion disc of a Be star, captured the same way at the closest approach of an eccentric orbit
  (the Be/X-ray binary).
* `overflow`: a donor that fills its lobe gives up its mass over `tuning.OVERFLOW_TIMESCALE_YR` (the low-mass X-ray
  binary).

Masses are solar masses, radii AU, axes AU, luminosities watts.
"""

import math

from planetgen import tuning
from planetgen.physics import constants
from planetgen.physics import unit_registry as units

_msun_to_kg = units.converter("Msun", "kg")

_YEAR_S = units.converter("yr", "s")(1.0)

PROTON_MASS_KG = 1.67262192369e-27  # [R4361]
THOMSON_CROSS_SECTION_M2 = 6.6524587321e-29  # [R4361]

MODES = ("wind", "disc", "overflow")
"""tuple: The ways the gas arrives."""

SPEED_OF_LIGHT_M_S = constants.SPEED_OF_LIGHT_M_S

_au_to_m = units.converter("au", "m")


def eddington_luminosity_w(mass_solar):
    """The Eddington luminosity of an accretor of `mass_solar` for hydrogen gas, `4 pi G M m_p c / sigma_T` (1.26e31 W per
    solar mass)."""
    return (4.0 * math.pi * constants.G * mass_solar * constants.SOLAR_MASS_TO_KG * PROTON_MASS_KG * SPEED_OF_LIGHT_M_S
            / THOMSON_CROSS_SECTION_M2)


def efficiency(kind, mass_solar, radius_au):
    """The share of accreted rest mass that is radiated: a black hole's `tuning.ACCRETION_EFFICIENCY_BLACK_HOLE`, a neutron
    star's surface potential `G M / (R c^2)` (about 0.17 for 1.4 solar masses on 12 km)."""
    if kind == "black-hole":
        return tuning.ACCRETION_EFFICIENCY_BLACK_HOLE
    return (constants.G * mass_solar * constants.SOLAR_MASS_TO_KG
            / (_au_to_m(radius_au) * SPEED_OF_LIGHT_M_S ** 2))


def donor_luminosity_lsun(mass_solar):
    """A hot massive donor's luminosity from its mass, `1.5 M^3.5` solar luminosities (a stand-in for the wind law's input;
    the member does not carry one)."""
    return 1.5 * mass_solar ** 3.5


def wind_mdot_kg_s(mass_solar):
    """The mass-loss rate of a hot massive star of `mass_solar`, kg/s."""
    rate_msun_yr = tuning.WIND_MASS_LOSS_MSUN_YR_AT_1E5_LSUN * (
        donor_luminosity_lsun(mass_solar) / 1.0e5) ** tuning.WIND_MASS_LOSS_LUMINOSITY_EXPONENT
    return _msun_to_kg(rate_msun_yr) / _YEAR_S


def wind_speed_m_s(mass_solar, radius_au):
    """The terminal wind speed: `tuning.WIND_SPEED_OVER_ESCAPE` times the surface escape speed."""
    escape = math.sqrt(2.0 * constants.G * mass_solar * constants.SOLAR_MASS_TO_KG / _au_to_m(radius_au))
    return tuning.WIND_SPEED_OVER_ESCAPE * escape


def captured_mdot_kg_s(accretor_mass_solar, axis_au, supply_kg_s, outflow_m_s, total_mass_solar):
    """The Bondi-Hoyle rate at which an accretor captures a flow of `supply_kg_s` leaving its source at `outflow_m_s`, at the
    separation `axis_au`, the pair moving at the circular speed of `total_mass_solar`."""
    a_m = _au_to_m(axis_au)
    v_orbit = math.sqrt(constants.G * total_mass_solar * constants.SOLAR_MASS_TO_KG / a_m)
    v_rel = math.sqrt(outflow_m_s ** 2 + v_orbit ** 2)
    gm = constants.G * accretor_mass_solar * constants.SOLAR_MASS_TO_KG
    return min(gm ** 2 * supply_kg_s / (a_m ** 2 * outflow_m_s * v_rel ** 3), supply_kg_s)


def overflow_mdot_kg_s(donor_mass_solar, donor_kind):
    """The rate a lobe-filling donor gives up mass, kg/s: its mass over `tuning.OVERFLOW_TIMESCALE_YR` of its kind (a
    main-sequence donor, else a giant's)."""
    timescale = tuning.OVERFLOW_TIMESCALE_YR["main-sequence" if donor_kind == "main-sequence" else "giant"]
    return _msun_to_kg(donor_mass_solar) / (timescale * _YEAR_S)


def x_ray_luminosity_w(mode, accretor_kind, accretor_mass_solar, accretor_radius_au, donor_kind, donor_mass_solar,
                       donor_radius_au, axis_au, eccentricity):
    """
    The X-ray luminosity in watts of a neutron star or a black hole that accretes in `mode` from a donor on an orbit of
    semi-major axis `axis_au` and `eccentricity`: `eta Mdot c^2`, never above the accretor's Eddington limit. The wind
    mode reads the mean separation, the disc mode the closest approach (the outburst).
    """
    if mode not in MODES:
        raise ValueError(f"unknown accretion mode {mode!r}; choose from {MODES}")
    total = accretor_mass_solar + donor_mass_solar
    if mode == "wind":
        mdot = captured_mdot_kg_s(accretor_mass_solar, axis_au, wind_mdot_kg_s(donor_mass_solar),
                                  wind_speed_m_s(donor_mass_solar, donor_radius_au), total)
    elif mode == "disc":
        mdot = captured_mdot_kg_s(accretor_mass_solar, axis_au * (1.0 - eccentricity),
                                  _msun_to_kg(tuning.BE_DISC_MASS_LOSS_MSUN_YR) / _YEAR_S,
                                  tuning.BE_DISC_OUTFLOW_KMS * 1000.0, total)
    else:
        mdot = overflow_mdot_kg_s(donor_mass_solar, donor_kind)
    power = efficiency(accretor_kind, accretor_mass_solar, accretor_radius_au) * mdot * SPEED_OF_LIGHT_M_S ** 2
    return min(power, eddington_luminosity_w(accretor_mass_solar))
