# planetgen/physics/formation.py

"""
Formation
=========

Protoplanetary disk physics: the snow line, the minimum-mass solar
nebula's surface density, and the isolation mass a planet can grow to.
"""

import math

from planetgen import tuning
from planetgen.physics import constants as physical_constants
from planetgen.util.checks import finite_domain


@finite_domain()
def snow_line_au(luminosity_w):
    """
    The snow line (ice condensation point, ~170K) for a star of a given
    luminosity, in AU -- real protoplanetary disk temperature falls off
    with stellar flux, i.e. with distance^-2, the same physical reasoning
    `calculate_habitable_zone` already uses for its own sqrt(luminosity)
    boundaries (see `physical_constants.SNOW_LINE_AU_AT_1_LSUN`).

    Args:
        luminosity_w (float): The star's luminosity, in Watts.

    Returns:
        float: The snow line's distance from the star, in AU.
    """
    solar_lum = luminosity_w / physical_constants.SOLAR_LUMINOSITY
    return physical_constants.SNOW_LINE_AU_AT_1_LSUN * math.sqrt(solar_lum)


@finite_domain()
def disk_surface_density_scale(star_mass_kg):
    """
    How much a star's own protoplanetary disk's solid surface density
    should be scaled relative to the Sun's (MMSN) value, based on the
    real, observed disk-mass-vs-stellar-mass relation -- see
    `tuning.DISK_MASS_STELLAR_MASS_EXPONENT`'s docstring for
    the literature basis. Feeds `mmsn_surface_density_gcm2`'s own
    `density_scale` argument.

    Args:
        star_mass_kg (float): The star's mass, in kilograms.

    Returns:
        float: A dimensionless scale factor, 1.0 for a solar-mass star.
    """
    solar_masses = star_mass_kg / physical_constants.SOLAR_MASS_TO_KG
    return solar_masses ** tuning.DISK_MASS_STELLAR_MASS_EXPONENT


@finite_domain()
def mmsn_surface_density_gcm2(distance_au, snow_line_au, density_scale=1.0):
    """
    The Minimum Mass Solar Nebula's solid surface density at a given
    distance from the star (Hayashi 1981) -- see
    `physical_constants.MMSN_SOLID_SURFACE_DENSITY_SOL_GCM2`'s docstring
    for the model and its snow-line ice-boost jump.

    Args:
        distance_au (float): Distance from the star, in AU.
        snow_line_au (float): This star's own snow line (see
                              `snow_line_au`), in AU.
        density_scale (float, optional): A star-dependent scale factor
                                         (see `disk_surface_density_scale`)
                                         applied on top of the Sun's own
                                         MMSN normalization. Defaults to
                                         1.0 (a solar-mass star's disk).

    Returns:
        float: The solid surface density at this distance, in g/cm^2.
    """
    density = (
        physical_constants.MMSN_SOLID_SURFACE_DENSITY_SOL_GCM2
        * density_scale
        * distance_au ** physical_constants.MMSN_SURFACE_DENSITY_EXPONENT
    )
    if distance_au >= snow_line_au:
        density *= physical_constants.SNOW_LINE_ICE_BOOST_FACTOR
    return density


@finite_domain()
def isolation_mass_kg(distance_au, surface_density_gcm2, star_mass_kg):
    """
    The oligarchic-growth isolation mass (Lissauer 1993; Kokubo & Ida
    2000, 2002): the mass a growing embryo reaches once it has cleared its
    own feeding zone of width `b` mutual Hill radii
    (`tuning.MUTUAL_HILL_RADII_SEPARATION`) -- the same `b`
    `StarSystem._mutual_min_distance_au` uses for adjacent-planet spacing,
    so this count estimate and that spacing rule are provably consistent.

    Derivation: M_iso = 2*pi*a*(b*R_H)*Sigma, where R_H = a*(M_iso /
    (3*M_star))^(1/3) is the embryo's *own* Hill radius (it hasn't met a
    neighbor yet, so this is the single-body form, not the mutual one).
    Substituting and solving for M_iso (it appears on both sides) gives
    the closed form below:

        M_iso = (2*pi * b * Sigma * a^2)^(3/2) / (3 * M_star)^(1/2)

    Args:
        distance_au (float): Distance from the star, in AU.
        surface_density_gcm2 (float): Local solid surface density at this
                                      distance (see
                                      `mmsn_surface_density_gcm2`), in
                                      g/cm^2.
        star_mass_kg (float): The star's mass, in kilograms.

    Returns:
        float: The isolation mass, in kilograms.
    """
    distance_m = distance_au * physical_constants.AU_TO_M
    surface_density_kgm2 = surface_density_gcm2 * 10  # 1 g/cm^2 = 10 kg/m^2
    b = tuning.MUTUAL_HILL_RADII_SEPARATION
    base = 2 * math.pi * b * surface_density_kgm2 * distance_m ** 2
    return base ** (3 / 2) / (3 * star_mass_kg) ** 0.5
