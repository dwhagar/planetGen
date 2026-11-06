# planetgen/galaxy/gravity_settings.py

"""
The numbers of the gravity map's hierarchy that its cache stamp reads (MAP.385, API.24; docs/design/orbits/gravity-map.md)
===========================================================================================================================

`galaxy.gravity_tree` (numba) builds the hierarchy and chooses the opening angle; the web app has to know whether a sector's cached grid is current without
importing it, and the stamp (`db.gravity_stamp.source_stamp`) names these numbers, so they live here, where both read them.
"""

LEAF_PC = 0.5
"""float: The edge of a leaf cell, parsecs: the finest grain the far sources are aggregated at (PERF.72: with whole 4 pc sectors as leaves a dense sector costs ten times more)."""
MIN_THETA = 0.3
"""float: The smallest opening angle `choose_theta` tries."""
MAX_THETA = 0.7
"""float: The largest opening angle that is allowed (1.0 gives 4 to 5 percent errors, PERF.72)."""
THETA_STEP = 0.05
"""float: The step `choose_theta` raises the opening angle by."""
PAIR_BUDGET = 50_000_000
"""int: The point-source pairs one target's list may cost with its points (about 0.1 s on four cores, PERF.72)."""
