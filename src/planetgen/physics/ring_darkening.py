# planetgen/physics/ring_darkening.py

"""
How bright a ring is, and the shade it throws (GEN.1055; docs/design/systems/weird-worlds.md section 9.2, rows R6 and R8)
========================================================================================================================

Two pure functions of the rings family (`rings.py`). Nothing here draws a number from a stream and nothing is stored: the ring painter (MAP.275) and the climate term of
the weird worlds (W1) read the results once a lane wires them.

R6, bright and dark rings. A ring is made of ice and of dark material (rock, dust, organics). Saturn's rings are 90 to 99 percent water ice and bright (Cuzzi et al. 2018, recalled),
Uranus's and Neptune's are dark; the Cassini ring mass and the flux of non-icy dust gave Saturn's rings an age of order 100 million years (Iess et al. 2019; the
age is contested, 2ND-HAND: no paper was read). Two things set the ice share a ring has now:

* **What it began with.** `rings.RANGES["ice_fraction"]`: a broad ring beyond the snow line is 90 to 99 percent ice, a faint or narrow one little, and a ring inside the snow line none.
* **What has fallen on it since.** Interplanetary dust darkens whatever it lands on. The dust that falls on a ring is a mass flux `F` (kg per m2 per s) and the ring holds a mass `Sigma`
  per m2, so after `t` seconds the ring has gained `F t / Sigma` of its own mass in dark material (`pollution_gain`). A thick ring (Saturn's B ring, 500 kg/m2 and more) is polluted by well under a percent in the
  age of the Solar System, a thin one (a ring of 1 kg/m2 at Uranus's distance) doubles its mass in dark dust and goes dark. The ring's mass per area is the optical depth over the cross-section
  the particles show: `Sigma = 4/3 rho_m r tau` (`rings.optical_depth` inverted), with the density of the mixture (`rings.mixture_density_kg_m3`).

The flux itself (`pollution_flux`) is the Saturn value scaled to the ring's distance from the star by the falling interplanetary dust density, and to the planet by gravitational focusing
(the dust is pulled in, `1 + (v_esc / v_inf)^2`).

The ice share left is `ice0 / (1 + gain)` (the added mass is all dark), and the albedo follows from the dark share by an intimate-mixture law: a ring of clean ice is bright, and
a few percent of dark material takes the brightness down fast, as the reflectance of grain mixtures does (Clark and Lucey 1984; recalled). `albedo` is the half-saturation form
`A = A_dark + (A_ice - A_dark) / (1 + d / d_half)`, with `d` the dark share.

R8, the ring's shadow. The ring lies in the planet's equatorial plane. In a season where the Sun stands at declination `delta` (the latitude under it) a point of the planet's winter
hemisphere sees the Sun through the ring when the line to the Sun crosses the equatorial plane between the ring's inner and outer radius: `shadow_transmission` finds, for a latitude
and a declination, the share of the day's sunlight that comes through (the path through the ring is longer by `1 / sin(delta)`, so the light left is `exp(-tau / sin(delta))`). The
shaded latitudes, the daily light left there and the winter latitude it is worst at are the inputs of the W1 climate term; the moons' shade follows from the same line.

Use:

    flux = pollution_flux(distance_au, planet_mass_kg, ring_radius_km)             # kg per m2 per s
    sigma = surface_density_kg_m2(optical_depth, particle_radius_m, ice_fraction)
    brightness = ring_brightness(ice0, age_yr, sigma, flux)                        # Brightness
    light = shadow_transmission(latitude_deg, declination_deg, inner_radii, outer_radii, optical_depth)   # 0..1
"""

import math
from dataclasses import dataclass

from planetgen.physics import constants as physical_constants, rings
from planetgen.util.checks import finite_domain

CLEAN_ICE_ALBEDO = 0.9
"""float: The reflectance of a ring of clean water-ice grains (the geometric albedo of fresh ice, 0.9 to 1 in the visible). [Stand-in RGD1, tier T6, tweak L3.]"""
DARK_ALBEDO = 0.04
"""float: The reflectance of a ring of dark material only (Uranus's and Neptune's rings, 0.03 to 0.05, recalled, 2ND-HAND). [Stand-in RGD2, tier T4, tweak L0.]"""
HALF_DARKENING_SHARE = 0.04
"""float: The dark share at which the ring is half way between clean ice and dark (an intimate mixture of ice and dark grains is half as bright at a few percent dark
material; Clark and Lucey 1984, recalled). Tuned so that Saturn's A and B rings (1 to 10 percent dark) come out at 0.3 to 0.8 and the C ring (about 30 percent) at 0.1. [Stand-in RGD1, tier T6, tweak L3.]"""

SATURN_POLLUTION_FLUX_KG_M2_S = 4.5e-17
"""float: The mass of non-icy micrometeoroids that falls on Saturn's rings per m2 and second, of order 1e-17 to 1e-16 (recalled, 2ND-HAND; the value that
gives Saturn's rings the non-icy 0.1 to 2 percent they show within a few hundred million years). [Stand-in RGD3, tier T3, tweak L0.]"""
REFERENCE_DISTANCE_AU = 9.58
"""float: Saturn's distance from the Sun, where `SATURN_POLLUTION_FLUX_KG_M2_S` holds."""
REFERENCE_PLANET_MASS_KG = 5.683e26
"""float: Saturn's mass, for the focusing of the reference flux."""
REFERENCE_RING_RADIUS_KM = 1.0e5
"""float: A radius in Saturn's rings (the B ring is at 92 000 to 118 000 km), for the focusing of the reference flux."""
FLUX_DISTANCE_EXPONENT = -1.3
"""float: How the interplanetary dust flux falls with distance from the star: the density of the dust falls as r^-1.3 (recalled; beyond Saturn it is
poorly measured). [Stand-in RGD4, tier T6, tweak L0.]"""
DUST_SPEED_AT_INFINITY_KMS = 5.0
"""float: The speed of the dust before the planet's pull is added (a few km/s for the dust of the outer Solar System). [Stand-in RGD4, tier T6, tweak L0.]"""

_G = physical_constants.G
_YEAR_S = physical_constants.SECONDS_PER_YEAR


def _focus(planet_mass_kg, ring_radius_km):
    """The gravitational focusing factor `1 + (v_esc / v_inf)^2` at the ring."""
    escape_kms = math.sqrt(2.0 * _G * planet_mass_kg / (ring_radius_km * 1.0e3)) / 1.0e3
    return 1.0 + (escape_kms / DUST_SPEED_AT_INFINITY_KMS) ** 2


@finite_domain()
def pollution_flux(distance_au, planet_mass_kg, ring_radius_km):
    """
    The mass of dark dust that falls on a ring per m2 and second, kg.

    Args:
        distance_au (float): The planet's distance from its star, AU.
        planet_mass_kg (float): The planet's mass, kg.
        ring_radius_km (float): A radius in the ring, km (the mean of its edges).
    """
    if distance_au <= 0 or planet_mass_kg <= 0 or ring_radius_km <= 0:
        raise ValueError("distance, mass and ring radius must be positive")
    reference = _focus(REFERENCE_PLANET_MASS_KG, REFERENCE_RING_RADIUS_KM)
    scale = (distance_au / REFERENCE_DISTANCE_AU) ** FLUX_DISTANCE_EXPONENT
    return SATURN_POLLUTION_FLUX_KG_M2_S * scale * _focus(planet_mass_kg, ring_radius_km) / reference


@finite_domain()
def surface_density_kg_m2(optical_depth, particle_radius_m, ice_fraction):
    """
    The ring's mass per m2 of the area it covers: `4/3 rho_m r tau` (the inverse of `rings.optical_depth`), kg.

    Args:
        optical_depth (float): The mean optical depth.
        particle_radius_m (float): The effective particle radius, m.
        ice_fraction (float): The ice share by mass (it sets the density of the mixture).
    """
    if optical_depth < 0 or particle_radius_m <= 0:
        raise ValueError("optical depth must not be negative and the particle radius must be positive")
    return 4.0 / 3.0 * rings.mixture_density_kg_m3(min(1.0, max(0.0, ice_fraction))) * particle_radius_m * optical_depth


@finite_domain()
def pollution_gain(age_yr, surface_density, flux_kg_m2_s):
    """
    The mass of dark dust a ring has gained since it began, as a share of its own mass: `F t / Sigma`.

    Args:
        age_yr (float): The time since the ring began, yr.
        surface_density (float): The ring's mass per m2, kg.
        flux_kg_m2_s (float): The dust flux, `pollution_flux`.
    """
    if age_yr < 0 or flux_kg_m2_s < 0:
        raise ValueError("age and flux must not be negative")
    if surface_density <= 0:
        raise ValueError("the ring's mass per area must be positive")
    return flux_kg_m2_s * age_yr * _YEAR_S / surface_density


@finite_domain()
def albedo(ice_fraction):
    """
    The reflectance of a ring that is `ice_fraction` water ice by mass and the rest dark: `A_dark + (A_ice - A_dark) / (1 + d / d_half)` with `d = 1 - ice_fraction`.
    """
    dark = 1.0 - min(1.0, max(0.0, ice_fraction))
    return DARK_ALBEDO + (CLEAN_ICE_ALBEDO - DARK_ALBEDO) / (1.0 + dark / HALF_DARKENING_SHARE)


@dataclass(frozen=True)
class Brightness:
    """
    How bright a ring is now.

    Attributes:
        ice_fraction (float): The ice share left after the dust that fell on it.
        dark_fraction (float): One minus `ice_fraction`.
        pollution_gain (float): The dark mass it gained, as a share of its own mass.
        albedo (float): Its reflectance, `albedo`.
    """

    ice_fraction: float
    dark_fraction: float
    pollution_gain: float
    albedo: float


def ring_brightness(ice_fraction, age_yr, surface_density, flux_kg_m2_s):
    """
    The ice share and reflectance of a ring of initial `ice_fraction` after `age_yr` of dust falling on it (R6).

    Args:
        ice_fraction (float): The ice share by mass when the ring began (0 to 1).
        age_yr (float): The ring's age, yr (`ring_age_yr`).
        surface_density (float): Its mass per m2, kg (`surface_density_kg_m2`).
        flux_kg_m2_s (float): The dust flux on it (`pollution_flux`).

    Returns:
        Brightness
    """
    if not 0.0 <= ice_fraction <= 1.0:
        raise ValueError("the ice share must be between 0 and 1")
    gain = pollution_gain(age_yr, surface_density, flux_kg_m2_s)
    ice = ice_fraction / (1.0 + gain)
    return Brightness(ice_fraction=ice, dark_fraction=1.0 - ice, pollution_gain=gain, albedo=albedo(ice))


@finite_domain()
def shadow_transmission(latitude_deg, declination_deg, inner_radii, outer_radii, optical_depth, longitude_steps=180):
    """
    The share of a day's sunlight that reaches a latitude of the planet through its ring (R8), 1 where the ring casts no shade.

    The ring is a thin annulus in the equatorial plane from `inner_radii` to `outer_radii` planet radii. A point on the planet at `latitude_deg` sees the Sun, whose subsolar
    latitude is `declination_deg`, through the ring when the line to the Sun crosses the equatorial plane at a distance from the planet's centre between the two edges. The light
    that comes through is `exp(-tau / sin(delta))` (the slant path through a ring of optical depth `tau`). The mean is over the lit half of a turn of the planet, weighted by the
    cosine of the Sun's zenith angle (the daily mean insolation) .

    Args:
        latitude_deg (float): The point's latitude, degrees (south negative).
        declination_deg (float): The Sun's declination relative to the ring plane, degrees (the planet's obliquity at a solstice, 0 at an equinox).
        inner_radii, outer_radii (float): The ring's edges in planet radii.
        optical_depth (float): The ring's mean optical depth.
        longitude_steps (int): How many longitudes sample the turn.

    Returns:
        float: The daily mean sunlight left, as a share of what a ringless planet would have there; 1.0 where the Sun is below the horizon all day.
    """
    if inner_radii <= 0 or outer_radii < inner_radii or optical_depth < 0 or longitude_steps < 8:
        raise ValueError("the ring's edges must be positive and ordered, its optical depth not negative")
    phi = math.radians(latitude_deg)
    delta = math.radians(declination_deg)
    sin_delta = math.sin(delta)
    cos_delta = math.cos(delta)
    lit = 0.0
    kept = 0.0
    for k in range(longitude_steps):
        lam = 2.0 * math.pi * (k + 0.5) / longitude_steps    # longitude from the subsolar meridian
        px = math.cos(phi) * math.cos(lam)
        py = math.cos(phi) * math.sin(lam)
        pz = math.sin(phi)
        weight = px * cos_delta + pz * sin_delta           # cosine of the zenith angle of the Sun
        if weight <= 0.0:
            continue
        lit += weight
        through = 1.0
        if abs(sin_delta) > 1.0e-9 and pz * sin_delta < 0.0:
            # The ray from the point to the Sun meets the equatorial plane at t = -pz / sin(delta).
            t = -pz / sin_delta
            hx = px + t * cos_delta
            hy = py
            if inner_radii <= math.hypot(hx, hy) <= outer_radii:
                through = math.exp(-optical_depth / abs(sin_delta))
        kept += weight * through
    return 1.0 if lit <= 0.0 else kept / lit
