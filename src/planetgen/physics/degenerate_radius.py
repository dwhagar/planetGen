"""
The radius of a cold degenerate star (a white dwarf, or the helium donor of an AM CVn binary) from its mass alone: the
zero-temperature mass-radius relation of Eggleton, as Verbunt and Rappaport (1988, ApJ 332, 193, eq. 15) and Rappaport,
Verbunt and Joss (1983) quote it [second-hand; the fit itself is unpublished]:

    R / R_sun = 0.0114 [(M / M_Ch)^(-2/3) - (M / M_Ch)^(2/3)]^(1/2) [1 + 3.5 (M / M_p)^(-2/3) + (M / M_p)^(-1)]^(-2/3)

with `M_Ch` 1.44 solar masses (the Chandrasekhar mass) and `M_p` 0.00057 solar masses (where the electrons stop being
degenerate). It holds from planetary masses to the Chandrasekhar mass and is the standard first reading for the
semi-degenerate donors of AM CVn systems, which it undersizes by up to a factor of two when the donor is still hot (Deloye
et al. 2007); the generator takes it as the donor's radius (GEN.1686). A carbon-oxygen white dwarf with a known
temperature uses `datasets.white_dwarf_models.radius_solar`, the cooling tracks, instead. Confident Unverified, L0.
"""

CHANDRASEKHAR_MASS_SOL = 1.44
"""float: `M_Ch`, solar masses, where the radius of the relation reaches zero."""

DEGENERACY_MASS_SOL = 0.00057
"""float: `M_p`, solar masses, below which the body is a planet held up by Coulomb pressure rather than degenerate electrons."""

RADIUS_SCALE_SOL = 0.0114
"""float: The leading coefficient, solar radii."""


def zero_temperature_radius_sol(mass_sol):
    """
    The radius, solar radii, of a cold degenerate body of `mass_sol` solar masses (see the module doc).

    Raises:
        ValueError: Unless `0 < mass_sol < CHANDRASEKHAR_MASS_SOL`.
    """
    if not 0.0 < mass_sol < CHANDRASEKHAR_MASS_SOL:
        raise ValueError(f"a degenerate body weighs between 0 and {CHANDRASEKHAR_MASS_SOL} solar masses, not {mass_sol}")
    ratio = mass_sol / CHANDRASEKHAR_MASS_SOL
    planet = mass_sol / DEGENERACY_MASS_SOL
    return (RADIUS_SCALE_SOL * (ratio ** (-2.0 / 3.0) - ratio ** (2.0 / 3.0)) ** 0.5
            * (1.0 + 3.5 * planet ** (-2.0 / 3.0) + planet ** -1.0) ** (-2.0 / 3.0))
