# planetgen/physics/main_sequence.py

"""
The main sequence of a star from one table (GEN.652): the MIST isochrones (`datasets/isochrones.py`), read at the
star's initial mass, its age and its birth [Fe/H].

A main-sequence star's luminosity, effective temperature and radius are one thing in the generator, whoever asks: the
star's own record (`stellar_evolution.evolve_star`), the pre-main-sequence tracks that contract onto it, the sampling
tables that need a star's brightness at the end of its main sequence, and the composition note
(`stellar_composition.star_makeup`) that reports the same star's table row. They all read `state`, `zams_state` and
`end_state` here. The mass-only piecewise laws they used to read (luminosity by mass in pieces, a radius by mass to a
power) are gone.

The phase times (`ms_lifetime`, Hurley, Pols and Tout 2000) are a second source that does not agree with MIST's own end
of the main sequence to the year: a star the lifetime keeps on the main sequence a little past MIST's edge, or not yet
on it, takes the row at that edge (`isochrones.main_sequence_point(..., clamp=True)`). A mass lighter than the grid's
lowest (0.1 solar masses; the generator makes stars from 0.08) is extended with the power laws of luminosity and radius
against mass at the grid's edge: a stand-in until a track reaches it.
"""

import functools
import math
from collections import namedtuple

from planetgen.datasets import isochrones
from planetgen import tuning

State = namedtuple("State", "luminosity_sol temperature_k radius_sol mass_sol")
"""A main-sequence star's luminosity (Lsun), effective temperature (K), radius (Rsun) and present mass (Msun: MIST's, which
has lost what the winds of a heavy star take from it on the main sequence, GEN.940)."""

_EDGE_STEP = 1.2
"""float: How far above the grid's edge mass the slope of the extension is read (a ratio of masses)."""


def _point(mass_sol, age_gy, fe_h):
    return isochrones.main_sequence_point(mass_sol, age_gy, fe_h, clamp=True)


def _edge_mass(mass_sol, age_gy, fe_h, direction):
    """The lightest (`direction` 1) or heaviest (-1) mass, a little inside the grid, that has a point."""
    low, high = tuning.MAIN_SEQUENCE_GRID_MASS_RANGE_SOL
    edge = low * 1.0001 if direction > 0 else high * 0.9999
    return edge if _point(edge, age_gy, fe_h) is not None else None


def _extended(mass_sol, age_gy, fe_h):
    """The point of a mass the grid does not hold: the power laws of luminosity and radius against mass read at the
    edge it is past, carried out to the mass (a stand-in, see the module note)."""
    low, high = tuning.MAIN_SEQUENCE_GRID_MASS_RANGE_SOL
    below = mass_sol < low
    edge = _edge_mass(mass_sol, age_gy, fe_h, 1 if below else -1)
    if edge is None:
        raise ValueError(f"the isochrone grid has no main-sequence point near {mass_sol} solar masses")
    inner = edge * _EDGE_STEP if below else edge / _EDGE_STEP
    at_edge, at_inner = _point(edge, age_gy, fe_h), _point(inner, age_gy, fe_h)
    span = math.log10(inner / edge)
    shift = math.log10(mass_sol / edge)
    log_l = at_edge[0] + (at_inner[0] - at_edge[0]) / span * shift
    log_r = at_edge[2] + (at_inner[2] - at_edge[2]) / span * shift
    log_teff = math.log10(tuning.SUN_EFFECTIVE_TEMPERATURE_K) + 0.25 * (log_l - 2.0 * log_r)
    return log_l, log_teff, log_r, mass_sol * at_edge[3] / edge       # a star under the grid has lost nothing worth counting


def state(mass_sol, age_gy, fe_h=0.0):
    """
    The luminosity, effective temperature and radius of a main-sequence star of initial mass `mass_sol` (solar masses),
    `age_gy` Gy old, of birth [Fe/H] `fe_h`: the MIST isochrone row, held to the mass's own span on the main sequence
    (a younger or older age takes the row at the edge). Raises `planetgen.datasets.fetch.DataMissing` when the
    isochrones were not downloaded.
    """
    point = _point(mass_sol, age_gy, fe_h)
    if point is None:
        point = _extended(mass_sol, age_gy, fe_h)
    log_l, log_teff, log_r, present_mass = point
    return State(10.0 ** log_l, 10.0 ** log_teff, 10.0 ** log_r, present_mass)


_EDGE_MEMO = {"token": None, "rows": {}}
"""dict: The rows `_edge` worked out, for the grid `token` names (the sampling tables and every star's record ask for the
same few masses at the same compositions a million times)."""

_EDGE_MEMO_LIMIT = 1 << 20
"""int: How many rows `_EDGE_MEMO` holds before it starts over."""


def _edge(mass_sol, fe_h, departure):
    token = isochrones.grid_token()
    memo = _EDGE_MEMO
    if memo["token"] != token or len(memo["rows"]) >= _EDGE_MEMO_LIMIT:
        memo["token"], memo["rows"] = token, {}
    key = (mass_sol, fe_h, departure)
    found = memo["rows"].get(key)
    if found is None:
        row = isochrones.edge_point(mass_sol, fe_h, departure)
        if row is None:
            found = state(mass_sol, 1.0e3 if departure else 0.0, fe_h)
        else:
            found = State(10.0 ** row[0], 10.0 ** row[1], 10.0 ** row[2], row[3])
        memo["rows"][key] = found
    return found


def zams_state(mass_sol, fe_h=0.0):
    """The star's state as it arrives on the main sequence: the first row its mass has (read from `isochrones.edge_point`,
    a table of the rows, for the many callers that ask)."""
    return _edge(mass_sol, fe_h, False)


def end_state(mass_sol, fe_h=0.0):
    """The star's state as it leaves the main sequence: the last row its mass has, carried on to the age it leaves (for
    a mass that never leaves in the grid's 20 Gy, its row at the oldest isochrone). A subgiant begins at this brightness
    and a supergiant brightens from it."""
    return _edge(mass_sol, fe_h, True)
