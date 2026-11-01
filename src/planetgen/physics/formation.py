# planetgen/physics/formation.py

"""
Formation
=========

Protoplanetary disk physics: the snow line, the minimum-mass solar
nebula's surface density, and the isolation mass a planet can grow to.
"""

import math

from planetgen import tuning
from planetgen.physics import condensate_shares, constants as physical_constants, unit_registry as units
from planetgen.util import draw
from planetgen.util.checks import finite_domain

_kg_to_msun = units.converter("kg", "Msun")


_au_to_m = units.converter("au", "m")


@finite_domain()
def snow_line_au(luminosity_w):
    """
    The snow line's blackbody position (where the starlight alone warms a
    blackbody to the temperature water ice condenses at, 182 K, [R2726: Table 10]: `condensate_shares.ICE_CONDENSATION_K`), in AU -- real protoplanetary disk temperature falls off with
    stellar flux, i.e. with distance^-2, the same physical reasoning
    `calculate_habitable_zone` already uses for its own sqrt(luminosity)
    boundaries. This is the centre of the scatter a star draws its own
    line from (`draw_snow_line_au`), not the star's line.

    Args:
        luminosity_w (float): The star's luminosity, in Watts.

    Returns:
        float: The blackbody snow line's distance from the star, in AU.
    """
    return condensate_shares.ice_line_au(luminosity_w, "H2O")


def draw_snow_line_au(luminosity_w, stream=None):
    """
    One star's own snow line, in AU (GEN.458): the blackbody position
    (`snow_line_au`) times 10 to a normal draw of `tuning.SNOW_LINE_SCATTER_DEX`
    dex, redrawn while it lies beyond `tuning.SNOW_LINE_SCATTER_MAX_SIGMA`
    sigma (a factor of 4 either way), which keeps it inside the factor 2 to 3
    that real lines move with the accretion history and the dust (Mulders et
    al. 2015). Drawn once for a star and stored with it.

    Args:
        luminosity_w (float): The star's luminosity, in Watts.
        stream (draw.Stream or None): The star's own stream to draw from;
            `None` draws from the stream bound now.

    Returns:
        float: The star's snow line, in AU.
    """
    centre = snow_line_au(luminosity_w)
    scatter = tuning.SNOW_LINE_SCATTER_DEX
    limit = tuning.SNOW_LINE_SCATTER_MAX_SIGMA
    source = stream or draw
    while True:
        z = source.gauss(0.0, 1.0)
        if abs(z) <= limit:
            return centre * 10.0 ** (scatter * z)
        draw.redrawn("snow_line.scatter")


@finite_domain()
def disk_surface_density_scale(star_mass_kg, metal_ratio=1.0):
    """
    How much a star's own protoplanetary disk's solid surface density
    should be scaled relative to the Sun's (MMSN) value, based on the
    real, observed disk-mass-vs-stellar-mass relation -- see
    `tuning.DISK_MASS_STELLAR_MASS_EXPONENT`'s docstring for
    the literature basis. Feeds `mmsn_surface_density_gcm2`'s own
    `density_scale` argument.

    Args:
        star_mass_kg (float): The star's mass, in kilograms.
        metal_ratio (float): The star's metal mass fraction over the Sun-like reference's, Z / Z0 (GEN.464,
            `bulk_draw.metal_ratio_to_solar`): the disk's solids scale with it. It is `10^[Fe/H]` for a thin-disc star
            of solar abundance ratios and 1.7 to 3 times that for a metal-poor high-alpha star, whose oxygen and
            magnesium are enhanced. 1 (solar) leaves the scale as it was.

    Returns:
        float: A dimensionless scale factor, 1.0 for a solar-mass star of
            solar metallicity.
    """
    solar_masses = _kg_to_msun(star_mass_kg)
    return solar_masses ** tuning.DISK_MASS_STELLAR_MASS_EXPONENT * metal_ratio


@finite_domain()
def mmsn_surface_density_gcm2(distance_au, snow_line_au, density_scale=1.0):
    """
    The Minimum Mass Solar Nebula's solid surface density at a given
    distance from the star [R1854] -- see
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
        density *= condensate_shares.WATER_ICE_SOLID_BOOST
    return density


@finite_domain()
def isolation_mass_kg(distance_au, surface_density_gcm2, star_mass_kg, separation=None):
    """
    The oligarchic-growth isolation mass ([R2709]; Kokubo & Ida
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
        separation (float, optional): The feeding zone's width `b` in mutual Hill radii; `tuning.MUTUAL_HILL_RADII_SEPARATION`
                                      when omitted.

    Returns:
        float: The isolation mass, in kilograms.
    """
    distance_m = _au_to_m(distance_au)
    surface_density_kgm2 = surface_density_gcm2 * 10  # 1 g/cm^2 = 10 kg/m^2
    b = tuning.MUTUAL_HILL_RADII_SEPARATION if separation is None else separation
    base = 2 * math.pi * b * surface_density_kgm2 * distance_m ** 2
    return base ** (3 / 2) / (3 * star_mass_kg) ** 0.5
