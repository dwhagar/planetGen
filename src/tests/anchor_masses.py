# tests/anchor_masses.py

"""
The mass of a scattered stellar-mass black hole in the anchor law's masses (GEN.554). The scatter's mass cut is in the anchor law (5 to 20 solar masses): each
hole keeps its rank in the law, and a place whose metals are not the disk's (`remnant_density.stretch_array`) stretches the span above 5 solar masses, so a hole of
rank 14 sits at 13.96 or 14.63 solar masses there. A test of "above the cut" or "under the cut" for these holes asks the anchor mass, which this gives back.
"""

import numpy as np

from planetgen import tuning
from planetgen.galaxy import remnant_density
from planetgen.generation import phenomenon_scatter as scatter

_COLUMNS = scatter.PHENOMENON_SCATTER_COLUMNS
_MASS, _KIND, _SUBTYPE = (_COLUMNS.index(name) for name in ("mass_solar", "kind", "subtype"))
_POSITION = tuple(_COLUMNS.index(f"position_{axis}_mpc") for axis in "xyz")


def anchor_mass(row, shape, center_pc=None):
    """The mass of a scatter row in the anchor law: a stellar-mass black hole's is its mass taken back from the span of its place, any other row's is its own.
    A ring sector's holes (`phenomenon_scatter.ring_remnant_rows`) take the span of the sector's center, `center_pc`; every other row's is the span at its own position."""
    mass = row[_MASS]
    if mass is None or (row[_KIND], row[_SUBTYPE]) != ("black-hole", "stellar"):
        return mass
    if center_pc is not None:
        x, y, z = (np.array([float(c)]) for c in center_pc)
    else:
        x, y, z = (np.array([row[column] / 1000.0]) for column in _POSITION)
    span = float(remnant_density.stretch_array(x, y, z, shape)[0])
    low = tuning.BLACK_HOLE_MASS_RANGE_SOLAR[0]
    return low + (mass - low) / span
