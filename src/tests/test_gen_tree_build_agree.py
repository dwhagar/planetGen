"""
GEN.1017: the tree of stars a system draws and the stars its build makes are the same stars. `StarSystem._drawn_plan` hands the primary's birth [Fe/H] to the
tree (`multistar.Placed.fe_h`, GEN.705) and `_add_companions` makes each companion from its member record with the same seed, birth mass, age and [Fe/H] (GEN.1009),
then `stored_trees.matches` holds the built star to the member's mass and radius. A mismatch leaves the plan unspent, so over a sample of systems that draw their own
companions there must be none.
"""

from planetgen.generation import stored_trees
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.util import draw

TYPES = ("A1V", "F5V", "G2V", "K3V", "M3V")
SEEDS_PER_TYPE = 8
"""int: Systems of each type: 40 in all, of which about 17 draw one or more companions (326 of 750 on the 11-28 census, 495 member checks, none off)."""


def test_the_tree_and_the_build_make_the_same_companions(monkeypatch):
    checks = []
    original = stored_trees.matches

    def watched(star, member):
        agreed = original(star, member)
        checks.append(agreed)
        return agreed

    monkeypatch.setattr(stored_trees, "matches", watched)
    for star_type in TYPES:
        for seed in range(3000, 3000 + SEEDS_PER_TYPE):
            draw.set_run_seed(seed)
            config = SystemConfig()
            config.STAR_TYPE = star_type
            config.PLANETS = False
            config.MOONS = False
            StarSystem(system_config=config)
    assert checks, "no system of the sample drew a companion, so nothing was compared"
    assert all(checks), f"{checks.count(False)} of {len(checks)} built companions are not the member the tree drew"
