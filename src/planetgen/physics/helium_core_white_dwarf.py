# planetgen/physics/helium_core_white_dwarf.py

"""
The Helium-Core White Dwarf of a Close Binary (GEN.996)
=======================================================

A star under about 2.3 solar masses climbs the red giant branch with a degenerate helium core, and its radius is set by the core
and hardly at all by the envelope on top of it (Refsdal and Weigert 1971; Webbink, Rappaport and Savonije 1983), Rappaport et
al. 1995's fit: `R = R0 Mc^4.5 / (1 + 4 Mc^4) + 0.5` solar radii for a core of `Mc` solar masses, with `R0` 5500 at Pop I
metallicity (Z 0.02) and 3300 at Pop II (Z 0.001). The core grows to about 0.45 solar masses before the helium flash. A giant in
a close pair fills its Roche lobe when its radius reaches the lobe's, and the envelope it then loses (stable overflow, or a
common envelope that is ejected) leaves the helium core it had at that moment: a helium-core white dwarf of 0.18 to 0.45 solar
masses, where the single-star relation (`remnant_mass.white_dwarf_mass_sol`) gives 0.55 or more. About a tenth of the white
dwarfs of the solar neighbourhood are such low-mass objects (Liebert et al. 2005, as quoted by Brown et al. 2011, ApJ 737, 23);
the 20 pc sample has 7.4 percent under 0.45 solar masses (Giammichele et al. 2012).

This module is the physics only: the radius of a giant with a given core, the core a lobe of a given size is reached at, and the
white dwarf a pair leaves. The pair's size and masses are the generator's own (the period law and the mass-ratio law), so the share
of white dwarfs that are helium-core comes out of them and is not a rate set here: with the period law of a solar-type star
(`multistar.PERIOD_LAWS_LOG10_DAYS`) the companions whose lobe is reached with a 0.2 to 0.45 solar mass core are about 14 percent
of them, and about 6 percent of the white dwarfs of 1 to 1.5 solar mass progenitors have one. Not wired yet: the tree records
a white dwarf's mass from its progenitor (`stored_trees.matches`), so the helium-core mass needs the member record of GEN.1006 P3.
"""

import math

from planetgen import tuning
from planetgen.physics import constants, roche

_LOG_Z_SUN = math.log10(0.02)


def giant_radius_rsun(core_mass_sol, fe_h=0.0):
    """The radius, solar radii, of a red giant with a helium core of `core_mass_sol` solar masses (Rappaport et al. 1995) at
    birth [Fe/H] `fe_h`. Rises with the core; 1.6 at 0.15, 24 at 0.3 and 154 at 0.45 solar masses for Pop I."""
    m = core_mass_sol
    return _radius_constant(fe_h) * m ** 4.5 / (1.0 + 4.0 * m ** 4) + 0.5


def _radius_constant(fe_h):
    """`R0` of the radius relation: straight in log Z between the two points of `tuning.GIANT_RADIUS_BY_CORE_MASS_RSUN`
    (Z 0.02 at [Fe/H] 0), held beyond them."""
    (low_logz, low), (high_logz, high) = tuning.GIANT_RADIUS_BY_CORE_MASS_RSUN
    logz = -(_LOG_Z_SUN + fe_h)
    return low + (high - low) * min(max((low_logz - logz) / (low_logz - high_logz), 0.0), 1.0)


def core_mass_at_radius_sol(radius_rsun, fe_h=0.0):
    """The helium core mass, solar masses, a giant has when its radius is `radius_rsun` solar radii, or None when that radius is
    reached before the core is `tuning.HELIUM_CORE_MIN_MASS_SOL` (a subgiant: not this relation) or only past the flash core
    `tuning.HELIUM_FLASH_CORE_MASS_SOL` (the star ignites helium first)."""
    low, high = tuning.HELIUM_CORE_MIN_MASS_SOL, tuning.HELIUM_FLASH_CORE_MASS_SOL
    if not giant_radius_rsun(low, fe_h) <= radius_rsun <= giant_radius_rsun(high, fe_h):
        return None
    for _ in range(60):
        middle = 0.5 * (low + high)
        if giant_radius_rsun(middle, fe_h) < radius_rsun:
            low = middle
        else:
            high = middle
    return 0.5 * (low + high)


def overflow_core_mass_sol(progenitor_sol, companion_sol, axis_au, eccentricity=0.0, fe_h=0.0):
    """
    The helium core, solar masses, of a giant of `progenitor_sol` solar masses at the moment it fills its Roche lobe in a pair
    with a star of `companion_sol` solar masses at semi-major axis `axis_au` and `eccentricity` (the lobe is filled at
    periastron, `axis * (1 - e)`), which is the mass of the helium white dwarf it leaves; None when the star does not leave one:
    a progenitor over `tuning.HELIUM_CORE_PROGENITOR_MAX_SOL` (its core is not degenerate), a pair too wide for the giant to
    reach its lobe before the helium flash, or one so close that the lobe is filled before the core is
    `tuning.HELIUM_CORE_MIN_MASS_SOL`.
    """
    if progenitor_sol > tuning.HELIUM_CORE_PROGENITOR_MAX_SOL:
        return None
    lobe_rsun = (axis_au * (1.0 - eccentricity) * roche.roche_lobe_fraction(progenitor_sol, companion_sol)
                 * constants.AU_M / constants.SOLAR_RADIUS_M)
    return core_mass_at_radius_sol(lobe_rsun, fe_h)


def overflow_axis_au(core_mass_sol, progenitor_sol, companion_sol, fe_h=0.0):
    """The circular axis, AU, at which a giant of `progenitor_sol` solar masses fills its lobe with a helium core of
    `core_mass_sol` solar masses beside a star of `companion_sol` solar masses (the inverse of `overflow_core_mass_sol` for
    e = 0)."""
    radius_au = giant_radius_rsun(core_mass_sol, fe_h) * constants.SOLAR_RADIUS_M / constants.AU_M
    return roche.separation_for_filling(radius_au, 1.0, progenitor_sol, companion_sol)
