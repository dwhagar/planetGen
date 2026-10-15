# planetgen/generation/limits.py

"""
Upper bounds on the generation inputs an admin can type: the Generate
pages (`planetgen/web/generate_page.py`, `planetgen/web/system_page.py`), the API
(`planetgen/api/routes.py`) and `planetgen`'s own argument parsing all check
against these same constants, so no form, request or command line can
ask for an unbounded run.

Each bound is far past anything a real run needs; it exists to turn an
absurd value (a typo, `9` * 25, a hostile request) into a clear error
instead of a job that runs for days or exhausts memory.
"""

from planetgen import tuning
from planetgen.galaxy import spherical
from planetgen.galaxy.skeleton import DEFAULT_MAX_SHELL
from planetgen.physics.units import pc_to_ly

MAX_GENERATE_RADIUS_PC = 200.0
"""float: The largest neighborhood radius (`--radius-pc`, the Generate
page's radius fields). A 200 pc sphere already holds about 520,000 sector
slots at the 4 pc sector standard; the default neighborhood (100 ly, about
31 pc) holds under 2,000."""

MAX_GENERATE_RADIUS_LY = pc_to_ly(MAX_GENERATE_RADIUS_PC)
"""float: `MAX_GENERATE_RADIUS_PC` in light years (about 652 ly), for
the API's generate-neighborhood `radius_ly`."""

MAX_GENERATE_SHELL = DEFAULT_MAX_SHELL
"""int: The highest shell index (`--shell`, `--max-shell`) -- the same cap
the density skeleton walks out to (real Milky-Way-scale galaxies end
around shell 3,800)."""

MAX_GENERATE_LIMIT = spherical.cells_in_band_range(MAX_GENERATE_SHELL, spherical.band_range(MAX_GENERATE_SHELL)[1])
"""int: The largest `--limit` for a shell batch: no shell up to
`MAX_GENERATE_SHELL` holds more cells than this (the cells grow with the
shell's surface, so the last one is the largest)."""

MAX_GENERATE_NEIGHBORHOODS = 100
"""int: The most random neighborhoods one run may ask for (`--neighborhoods`, GEN.97)."""

MAX_NUM_ORBITS = tuning.ABSOLUTE_MAX_SYSTEM_OBJECTS
"""int: The largest forced orbital slot count (`--num-orbits`, the API's
`num_orbits`) -- the generator's own ceiling on objects in a system."""
