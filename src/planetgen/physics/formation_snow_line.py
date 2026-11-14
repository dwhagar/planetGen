# planetgen/physics/formation_snow_line.py

"""
The snow line a star's disk was built round (GEN.1611)
=======================================================

A planet's ices and a belt's kind depend on whether it formed beyond the ice line of the star's disk, and the ice line moves with the star's luminosity
while the disk is there (a few Myr, `FORMATION_DISK_AGE_GY`). `formation.draw_snow_line_au` was given the star's luminosity today, so a white dwarf (0.0001 of a
Sun) got a snow line of 0.02 AU and a K giant 13 to 22 AU, and every "formed beyond the snow line" test round an evolved star read the wrong line.

`formation_luminosity_sol` is the luminosity the disk saw, the same for every star whatever it has become: that of its initial mass at `tuning.FORMATION_DISK_AGE_GY`
(`birth_luminosity.birth_luminosity_sol`: the BHAC15 isochrones to 1.4 solar masses, the zero-age main sequence above). A Sun-like star's disk saw 0.46 of its present
light and an M dwarf's several times its main-sequence light; an evolved star's disk saw its progenitor's, not the giant's or the white dwarf's. `formation_snow_line_au` is the
blackbody line of that luminosity and `draw_formation_snow_line_au` puts the star's own scatter on it (`formation.draw_snow_line_au`).

The water line (GEN.1678). The luminosity line is the same position at every stellar mass in units of the habitable zone, so a habitable-zone planet drew the same water class round a
0.2 solar-mass star as round the Sun; the formation papers find the opposite (Mulders 2015 [R3193], Raymond 2007 [R3683]: low-mass stars give dry habitable-zone planets), because a disk's ice line
is set by viscous heating and moves with the accretion rate, which falls steeply with the star's mass. `viscous_snow_line_au` is that line (Eq. 2 and 4 of [R3193], `tuning.VISCOUS_DISK_*`);
`water_snow_line_au` is the line the delivered water of a rocky body (`delivered_water`, the atmosphere recipes b1 to b4) is read against: the star's own drawn line (its scatter kept) moved by the
ratio of the viscous line to the luminosity line of its initial mass. Everything else (the planet chain, the belts, where a body is icy) still reads the luminosity line.
"""

from planetgen import tuning
from planetgen.physics import birth_luminosity, constants, formation, mass_loss
from planetgen.util import detmath

LUPUS_LOWEST_MASS_SOL = 0.03
"""float: The lowest stellar mass, in suns, the accretion-rate law was fitted to ([R3193] Eq. 3: 0.03 to 1.0); a lighter star is read at this mass."""


def formation_luminosity_sol(initial_mass_sol, fe_h):
    """The luminosity, in suns, of a star of this initial mass and [Fe/H] while its disk was there (`birth_luminosity.birth_luminosity_sol` at `tuning.FORMATION_DISK_AGE_GY`)."""
    return birth_luminosity.birth_luminosity_sol(initial_mass_sol, fe_h, tuning.FORMATION_DISK_AGE_GY)


def formation_snow_line_au(initial_mass_sol, fe_h):
    """The blackbody snow line, in AU, of `formation_luminosity_sol` (the star's own scatter, `draw_formation_snow_line_au`, goes on top of it)."""
    return formation.snow_line_au(formation_luminosity_sol(initial_mass_sol, fe_h) * constants.SOLAR_LUMINOSITY)


def draw_formation_snow_line_au(initial_mass_sol, fe_h, stream=None):
    """The star's own snow line, in AU, drawn once around `formation_snow_line_au` (`formation.draw_snow_line_au`: 0.2 dex, cut at 3 sigma)."""
    return formation.draw_snow_line_au(formation_luminosity_sol(initial_mass_sol, fe_h) * constants.SOLAR_LUMINOSITY, stream)


def initial_mass_sol_of(star):
    """The mass, in suns, `star` had when its disk formed: its stored initial mass, or where it stores none the mass its planets formed round (`mass_loss.member_formation_mass_kg`: a white dwarf's progenitor's)."""
    stored = getattr(star, "initial_mass_sol", None)
    if stored:
        return stored
    return mass_loss.member_formation_mass_kg(star) / constants.SOLAR_MASS_TO_KG


def draw_for_star(star, stream=None):
    """The snow line, in AU, a star's disk was built round (`draw_formation_snow_line_au` of its initial mass and its [Fe/H], 0 where it stores none)."""
    return draw_formation_snow_line_au(initial_mass_sol_of(star), getattr(star, "metallicity_fe_h", None) or 0.0, stream)


def viscous_snow_line_au(initial_mass_sol, opacity_cm2_g=None):
    """
    The viscous-disk snow line, in AU, of a star of this initial mass ([R3193] Eq. 2 with the accretion rate of Eq. 4 at `tuning.VISCOUS_DISK_AGE_MYR`):
    R = 2.1 au (M / Msun)^(1/3) (Mdot / 1e-8)^(4/9) (kappa / 770)^(2/9) (f / 100)^(-2/9) (alpha / 0.01)^(-2/9) (T_ice / 160 K)^(-10/9). Above 1 solar mass the accretion law is taken past
    where it was fitted (Lupus, 0.03 to 1.0), a stand-in extension. `opacity_cm2_g` is `tuning.VISCOUS_DISK_OPACITY_CM2_G` unless given.
    """
    opacity = tuning.VISCOUS_DISK_OPACITY_CM2_G if opacity_cm2_g is None else opacity_cm2_g
    mass = max(initial_mass_sol, LUPUS_LOWEST_MASS_SOL)
    mdot = (tuning.VISCOUS_DISK_MDOT_PREFACTOR_SOL_YR * mass ** tuning.VISCOUS_DISK_MDOT_MASS_EXPONENT
            * tuning.VISCOUS_DISK_AGE_MYR ** tuning.VISCOUS_DISK_MDOT_TIME_EXPONENT)
    return (2.1 * mass ** (1.0 / 3.0) * (mdot / 1.0e-8) ** (4.0 / 9.0) * (opacity / 770.0) ** (2.0 / 9.0)
            * (tuning.VISCOUS_DISK_GAS_TO_DUST / 100.0) ** (-2.0 / 9.0) * (tuning.VISCOUS_DISK_ALPHA / 0.01) ** (-2.0 / 9.0) * (tuning.VISCOUS_DISK_ICE_K / 160.0) ** (-10.0 / 9.0))


def _single_water_snow_line_au(star):
    stored = getattr(star, "snow_line_au", None)
    if not stored:
        return stored
    mass = initial_mass_sol_of(star)
    mean_au = formation_snow_line_au(mass, getattr(star, "metallicity_fe_h", None) or 0.0)
    return stored * viscous_snow_line_au(mass) / mean_au


def water_snow_line_au(star):
    """
    The snow line, in AU, the delivered water of a rocky body round `star` is read against (GEN.1678): the star's stored line times the viscous line over the luminosity line of its initial
    mass, so the star's own scatter stays. A pair adds its members' lines as the pair's stored line does (`detmath.hypot`). `None` for a star with no stored line.
    """
    members = getattr(star, "stars", None)
    if members and len(members) > 1:
        lines = [_single_water_snow_line_au(member) for member in members]
        return None if any(not line for line in lines) else detmath.hypot(*lines)
    return _single_water_snow_line_au(star)
