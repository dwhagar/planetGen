# planetgen/physics/formation_snow_line.py

"""
The snow line a star's disk was built round (GEN.1611, slice 0)
================================================================

A planet's ices and a belt's kind depend on whether it formed beyond the ice line of the star's disk, and the ice line moves with the star's luminosity
while the disk is there (a few Myr, `FORMATION_DISK_AGE_GY`). `formation.draw_snow_line_au` was given the star's luminosity today, so a white dwarf (0.0001 of a
Sun) got a snow line of 0.02 AU and a K giant 13 to 22 AU, and every "formed beyond the snow line" test round an evolved star read the wrong line.

`formation_luminosity_sol` is the luminosity the disk saw: the present one for a star still on the main sequence or contracting onto it (its line is what the
generator is calibrated on, and what the disk saw within a factor of a few), and for an evolved star (a subgiant, giant, supergiant, white dwarf) the luminosity its track had at
`FORMATION_DISK_AGE_GY`. `formation_snow_line_au` is the blackbody line of that luminosity. Pure: no draw (the scatter of `formation.draw_snow_line_au` goes on top of
it when slice 1 wires this), and not read by the generator yet (wiring register).
"""

from planetgen import tuning
from planetgen.physics import constants, formation, stellar_evolution
from planetgen.util import draw

UNEVOLVED_CLASSES = ("V", "VI", "PMS", "protostar")
"""tuple: The Yerkes classes of a star whose luminosity today is the one its disk saw (the main sequence, subdwarfs, and the contracting stages)."""


def formation_luminosity_sol(initial_mass_sol, fe_h, present_luminosity_sol, yerkes_class):
    """
    The luminosity, in suns, of the star while its disk was there.

    Args:
        initial_mass_sol (float): The star's mass when it formed (not a white dwarf's mass today).
        fe_h (float): Its [Fe/H].
        present_luminosity_sol (float): Its luminosity today.
        yerkes_class (str): Its Yerkes class today (`stellar_evolution.yerkes_class_of`).

    Returns:
        float: `present_luminosity_sol` for a class of `UNEVOLVED_CLASSES`; else the track's luminosity at `tuning.FORMATION_DISK_AGE_GY` (the present one if the track gives none).
    """
    if yerkes_class in UNEVOLVED_CLASSES:
        return present_luminosity_sol
    state = stellar_evolution.evolve_star(initial_mass_sol, tuning.FORMATION_DISK_AGE_GY, rng=draw.Stream(0), fe_h=fe_h)
    return present_luminosity_sol if state is None else state["luminosity_sol"]


def formation_snow_line_au(initial_mass_sol, fe_h, present_luminosity_sol, yerkes_class):
    """The blackbody snow line, in AU, of `formation_luminosity_sol` (the star's own scatter, `formation.draw_snow_line_au`, goes on top of it)."""
    return formation.snow_line_au(formation_luminosity_sol(initial_mass_sol, fe_h, present_luminosity_sol, yerkes_class) * constants.SOLAR_LUMINOSITY)
