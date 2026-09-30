# stellarObjects/generationLimits.py

"""
Upper bounds on the generation inputs an admin can type: the Generate
pages (`html/web/generate_page.py`, `html/web/system_page.py`), the API
(`html/api/routes.py`) and `generate.py`'s own argument parsing all check
against these same constants, so no form, request or command line can
ask for an unbounded run.

Each bound is far past anything a real run needs; it exists to turn an
absurd value (a typo, `9` * 25, a hostile request) into a clear error
instead of a job that runs for days or exhausts memory.
"""

from . import program_constants
from .galaxyGeometry import ring_sector_count
from .galaxySkeleton import DEFAULT_MAX_RING
from .utils import pc_to_ly

MAX_GENERATE_RADIUS_PC = 200.0
"""float: The largest neighborhood radius (`--radius-pc`, the Generate
page's radius fields). A 200 pc sphere already holds about 520,000 sector
slots at the 4 pc sector standard; the default neighborhood (100 ly, about
31 pc) holds under 2,000."""

MAX_GENERATE_RADIUS_LY = pc_to_ly(MAX_GENERATE_RADIUS_PC)
"""float: `MAX_GENERATE_RADIUS_PC` in light years (about 652 ly), for
the API's generate-neighborhood `radius_ly`."""

MAX_GENERATE_RING = DEFAULT_MAX_RING
"""int: The highest ring index (`--ring`, `--max-ring`) -- the same cap
the density skeleton walks out to (real Milky-Way-scale galaxies end
around ring 3,900)."""

MAX_GENERATE_LIMIT = ring_sector_count(MAX_GENERATE_RING)
"""int: The largest `--limit` for a ring batch: no ring up to
`MAX_GENERATE_RING` holds more slots than this."""

MAX_NUM_ORBITS = program_constants.ABSOLUTE_MAX_SYSTEM_OBJECTS
"""int: The largest forced orbital slot count (`--num-orbits`, the API's
`num_orbits`) -- the generator's own ceiling on objects in a system."""
