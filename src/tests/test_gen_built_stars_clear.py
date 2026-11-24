# tests/test_gen_built_stars_clear.py

"""
A system of three or more stars never holds one star inside another. A drawn tree is built from each member's own record (GEN.1009),
so the stars as built are the ones the orbits were drawn around; the independent verifier (`physics/clearance.py`) finds none inside
another or inside its Roche lobe, with no redraw of the companions to make it so (GEN.1015).
"""

import pytest

from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.util import draw
from tests.clearance_rows import hard_violations, rows_of_system


@pytest.mark.parametrize("seed", [241, 297, 59])
def test_a_forced_triple_the_first_draw_built_a_star_inside_another_clears(seed):
    config = SystemConfig()
    config.PLANETS = False
    config.STARS = 3
    with draw.counting() as redraws, draw.bound(seed):
        system = StarSystem(system_config=config)
    assert len(system.stars) == 3
    assert not hard_violations(rows_of_system(system))
    assert not [kind for kind in redraws if kind.endswith("overlap")], redraws
