# planetgen/galaxy/resolution.py

"""
How deep a cell is resolved (GEN.281)
=====================================

A cell that expects more systems than the **cell star ceiling** is filled only with the stars at or above the **fill
luminosity cap**; the rest of it is *unresolved* (docs/design/realistic-galactic-core.md section 4.2). Both are
settings of the galaxy (`store.GalaxySkeletonInfo.cell_star_ceiling` and `fill_luminosity_cap_sol`). Ordinary space
never reaches the ceiling, so it is untouched.

`cell_resolution` is a pure function of a cell's profile and those settings: the fill, the sector page and the API read
the same answer from the same address, so the unresolved count and mass are computed and never stored. It is a step and
not a slope: a cell expecting just under the ceiling is filled whole and one expecting just over it keeps only the
bright stars.

The ceiling is also the most a cell makes: when the stars at or above the fill luminosity cap would still be more than
the ceiling (the nuclear cluster's cell expects 15 million systems, 0.27% of them 40,000), the cap rises, along a
fixed ladder of `CAP_LADDER_STEP` steps so the bright-star tables are shared between cells, to the lowest rung whose
stars fit.
"""

import math
from dataclasses import dataclass

from planetgen import tuning
from planetgen.generation.star_population import placed_star_fraction


@dataclass(frozen=True)
class CellResolution:
    """
    What a cell resolves.

    Attributes:
        capped (bool): Whether the cell expects more systems than the ceiling, so only the bright ones are made.
        ceiling_systems (float): The cell star ceiling.
        luminosity_cap_sol (float): The fill luminosity cap, solar luminosities.
        expected_systems (float): The mean of the cell's Poisson draw, its volume allowed for.
        resolved_systems (float): The part of it made: all of it below the ceiling, else the stars at or above the cap
            (and those born at the scatter's mass limit or more).
        unresolved_systems (float): `expected_systems` less `resolved_systems`.
        unresolved_mass_solar (float): The mass of the unresolved systems, at `tuning.MEAN_SYSTEM_MASS_SOLAR` each.
    """
    capped: bool
    ceiling_systems: float
    luminosity_cap_sol: float
    expected_systems: float
    resolved_systems: float
    unresolved_systems: float
    unresolved_mass_solar: float


def expected_in_cell(profile):
    """The mean of a cell's Poisson draw: the profile's `expected_systems` (at an edge cube's volume) scaled to the
    cell's own volume, which differs from the cube's in the core ball and near the poles."""
    return profile.expected_systems * profile.volume_pc3 / profile.edge_pc ** 3


CAP_LADDER_STEP = 1.25
"""float: The ratio between neighbouring rungs of the ladder a raised fill luminosity cap climbs."""
CAP_LADDER_TOP_SOL = 1.0e9
"""float: The highest rung: no star of any population is this bright, so a cell always fits by here."""


def _resolved_at(profile, expected, cap_sol, mass_limit_sol):
    total = sum(profile.mix.values())
    if total <= 0.0:
        return 0.0
    share = sum(density * placed_star_fraction(cap_sol, mass_limit_sol, population)
                for population, density in profile.mix.items()) / total
    return expected * min(share, 1.0)


def cell_resolution(profile, ceiling_systems, luminosity_cap_sol, mass_limit_sol=None):
    """
    What a cell resolves under the ceiling and the cap. A cell over the ceiling is filled down to the fill luminosity
    cap, or to the lowest rung above it (`CAP_LADDER_STEP` apart) that makes no more than the ceiling.

    Args:
        profile (galaxy.profile.SectorProfile): The cell's profile.
        ceiling_systems (float): The cell star ceiling.
        luminosity_cap_sol (float): The fill luminosity cap (at least the brightest white dwarf, as every luminosity
            threshold is); the lowest the cell's own cap goes.
        mass_limit_sol (float, optional): The mass the galaxy's star scatter places every star from, `None` when it
            has no mass pass: those stars are resolved too, whatever their luminosity.

    Returns:
        CellResolution
    """
    expected = expected_in_cell(profile)
    if expected <= ceiling_systems:
        return CellResolution(False, float(ceiling_systems), float(luminosity_cap_sol), expected, expected, 0.0, 0.0)
    cap = float(luminosity_cap_sol)
    resolved = _resolved_at(profile, expected, cap, mass_limit_sol)
    if resolved > ceiling_systems:
        top = math.ceil(math.log(CAP_LADDER_TOP_SOL / cap, CAP_LADDER_STEP))
        low, high = 0, max(top, 1)
        while low < high:
            middle = (low + high) // 2
            if _resolved_at(profile, expected, cap * CAP_LADDER_STEP ** middle, mass_limit_sol) <= ceiling_systems:
                high = middle
            else:
                low = middle + 1
        cap = cap * CAP_LADDER_STEP ** low
        resolved = _resolved_at(profile, expected, cap, mass_limit_sol)
    unresolved = expected - resolved
    return CellResolution(True, float(ceiling_systems), cap, expected, resolved, unresolved,
                          unresolved * tuning.MEAN_SYSTEM_MASS_SOLAR)
