# planetgen/physics/birth_luminosity.py

"""
The luminosity of a star while its disk was there (GEN.553's birth step, for GEN.1611)
=======================================================================================

A planet's ices and a belt's kind depend on whether it formed beyond the snow line, and that line followed the luminosity the star had in its first few Myr, which is not the
luminosity it has today (an evolved star) and not its main-sequence luminosity either (a star still contracting onto the main sequence: a Sun-like star at 10 Myr is at about 0.46
of its present light, an M dwarf of 0.3 solar masses at 0.06 against 0.02 on the main sequence).

`birth_luminosity_sol(initial_mass_sol, fe_h, age_gy)` is that luminosity, in suns, at the star's age `age_gy` (the caller passes the age its disk is read at, GEN.1611's `FORMATION_DISK_AGE_GY`, 0.01 Gy; it is
not repeated here so that there is one stand-in for it). Pure: no draw.

* Where the BHAC15 isochrones hold the mass (0.01 to 1.4 solar masses, [R0254], `datasets.baraffe_isochrones`) the luminosity is theirs, interpolated in log mass and log age at
  that age (verified, solar metallicity, so `fe_h` does not enter). Below the lightest mass of the table the lightest row stands (the snow line of such an object is under the
  floor of the disk model anyway). Ages younger than the grid's first (0.5 Myr) read the first.
* Above 1.4 solar masses `stellar_evolution.evolve_star` gives the star's state at the age asked or just after its arrival on the main sequence, whichever is earlier: the
  zero-age luminosity of a star that has settled by then, so a star that has died by that age (30 solar masses live 6 Myr) is read as it began. The model track is brighter than BHAC15 at 10 Myr from 0.5 solar masses up (a Sun-like star: 0.75 against 0.46), so the two meet with a
  step at 1.4 solar masses (3.35 against 1.95 solar luminosities, a factor of 1.3 in the snow line) that this function does not hide; it is the model's step, listed on GEN.553.
"""

from planetgen.datasets import baraffe_isochrones
from planetgen.physics import pre_main_sequence, stellar_evolution
from planetgen.util import draw

ARRIVAL_MARGIN = 1.05
"""float: The factor on the contraction time at which a star above the table is read (just past its arrival on the main sequence). A numerical stand-in, tier T9, tweak L0: any margin that puts the read past the arrival gives the same zero-age light."""


def table_luminosity_sol(initial_mass_sol, age_gy):
    """The BHAC15 luminosity, in suns, of a star of this mass at this age, or None where the table does not hold the mass (above 1.4 solar masses, or the age is past 10 Gyr)."""
    grid = baraffe_isochrones.isochrones()
    mass = max(initial_mass_sol, float(grid[0].mass_solar[0]))   # the lightest row stands below the table
    age = max(age_gy, grid[0].age_gyr)
    state = baraffe_isochrones.interpolate(mass, age, grid)
    return None if state is None else 10.0 ** state["log_luminosity_solar"]


def birth_luminosity_sol(initial_mass_sol, fe_h, age_gy):
    """
    The luminosity, in suns, of a star of this initial mass and [Fe/H] at the age its disk saw (module docstring).

    Args:
        initial_mass_sol (float): The star's mass when it formed.
        fe_h (float): Its [Fe/H] (read above 1.4 solar masses only; BHAC15 is solar).
        age_gy (float): The age to read, Gy (the age the star's disk is read at).

    Returns:
        float: The luminosity, in suns.

    Raises:
        ValueError: The mass or the age is not above zero.
        fetch.DataMissing: The BHAC15 file is not downloaded.
    """
    if initial_mass_sol <= 0.0 or age_gy <= 0.0:
        raise ValueError(f"a birth luminosity needs a mass and an age above zero, not {initial_mass_sol} and {age_gy}")
    found = table_luminosity_sol(initial_mass_sol, age_gy)
    if found is not None:
        return found
    # Above the table a star is on the main sequence within 10 Myr (the model's contraction time is 8.6 Myr at 1.4 solar masses, 0.1 Myr at 10): its light is the
    # zero-age one, read just after its arrival, so a heavy star that has died or left the main sequence by the age asked (30 solar masses live 6 Myr) is read as it began.
    arrival_gy = pre_main_sequence.contraction_time_yr(initial_mass_sol, fe_h) * ARRIVAL_MARGIN * 1.0e-9
    state = stellar_evolution.evolve_star(initial_mass_sol, min(age_gy, arrival_gy), rng=draw.Stream(0), fe_h=fe_h)
    return state["luminosity_sol"]
