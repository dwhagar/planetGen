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
"""

from planetgen import tuning
from planetgen.physics import birth_luminosity, constants, formation, mass_loss


def formation_luminosity_sol(initial_mass_sol, fe_h):
    """The luminosity, in suns, of a star of this initial mass and [Fe/H] while its disk was there (`birth_luminosity.birth_luminosity_sol` at `tuning.FORMATION_DISK_AGE_GY`)."""
    return birth_luminosity.birth_luminosity_sol(initial_mass_sol, fe_h, tuning.FORMATION_DISK_AGE_GY)


def formation_snow_line_au(initial_mass_sol, fe_h):
    """The blackbody snow line, in AU, of `formation_luminosity_sol` (the star's own scatter, `draw_formation_snow_line_au`, goes on top of it)."""
    return formation.snow_line_au(formation_luminosity_sol(initial_mass_sol, fe_h) * constants.SOLAR_LUMINOSITY)


def draw_formation_snow_line_au(initial_mass_sol, fe_h, stream=None):
    """The star's own snow line, in AU, drawn once around `formation_snow_line_au` (`formation.draw_snow_line_au`: 0.2 dex, cut at the sourced spread of a factor of 3)."""
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
