"""GEN.1457: a Jupiter-family comet belongs to the star's controlling giant (`planetgen/physics/comet_control.py`, `Comet(giant_au=...)`).
The Solar System's 67P below is a test target only: the rule must admit it, and nothing in the module is centered on it."""

import math

import pytest

from planetgen import tuning
from planetgen.generation.comet import Comet
from planetgen.generation.config import SystemConfig
from planetgen.generation.redraw_kinds import KINDS
from planetgen.physics import comet_control
from planetgen.util import draw

LOW, HIGH = comet_control.JUPITER_FAMILY_TISSERAND


def _comet(seed, giant_au=None, **kwargs):
    draw.set_run_seed(seed)
    return Comet(SystemConfig(), primary_mass_solar=1.0, orbit_type="elliptical", giant_au=giant_au, **kwargs)


@pytest.fixture(autouse=True)
def _keep_draw_state():
    state = draw.getstate()
    yield
    draw.setstate(state)


def test_the_tisserand_parameter_of_a_known_orbit():
    # 67P: a = 3.46 AU, e = 0.64, i = 7.04 degrees about a giant at 5.2 AU is a Jupiter-family value.
    q = 3.46 * (1 - 0.64)
    assert comet_control.tisserand(5.2, q, 0.64, 7.04) == pytest.approx(2.75, abs=0.01)
    assert comet_control.controls(5.2, q, 0.64, 7.04)
    # A circular orbit in the giant's own plane at the giant's distance has T = 3 (the edge), a retrograde one is far below 2.
    assert comet_control.tisserand(5.2, 5.2, 0.0, 0.0) == pytest.approx(3.0)
    assert not comet_control.controls(5.2, q, 0.64, 170.0)


def test_the_controlling_giant_is_the_most_massive_one():
    class Body:
        def __init__(self, body_type, mass, distance):
            self.body_type, self.mass, self.distance = body_type, mass, distance

    planets = [Body('t', 9e27, 1.0), Body('g', 1e27, 2.0), Body('g', 3e27, 7.0), Body('g', 2e27, 11.0)]
    assert comet_control.controlling_giant_au(planets) == 7.0
    assert comet_control.controlling_giant_au(planets[:2][:1]) is None
    assert comet_control.controlling_giant_au([]) is None


def test_a_star_without_a_giant_draws_the_comets_it_drew_before():
    for seed in range(60):
        free = _comet(seed)
        none = _comet(seed, giant_au=None)
        assert (none.period_class, none.perihelion_distance_au, none.eccentricity, none.inclination_deg) == \
               (free.period_class, free.perihelion_distance_au, free.eccentricity, free.inclination_deg)


def test_a_giant_leaves_every_other_family_as_it_was():
    other = 0
    for seed in range(200):
        free = _comet(seed)
        held = _comet(seed, giant_au=5.2)
        if free.period_class != "jupiter_family":
            other += 1
            assert (held.period_class, held.perihelion_distance_au, held.eccentricity, held.inclination_deg) == \
                   (free.period_class, free.perihelion_distance_au, free.eccentricity, free.inclination_deg)
    assert other > 30


def test_a_jupiter_family_comet_about_a_giant_has_a_tisserand_value_in_the_range():
    seen = 0
    for seed in range(300):
        comet = _comet(seed, giant_au=5.2)
        if comet.period_class == "jupiter_family":
            seen += 1
            assert comet_control.controls(5.2, comet.perihelion_distance_au, comet.eccentricity, comet.inclination_deg)
    assert seen > 50


def test_a_giant_moves_the_perihelion_distances_of_the_family():
    free_q = [c.perihelion_distance_au for c in (_comet(s) for s in range(400)) if c.period_class == "jupiter_family"]
    held_q = [c.perihelion_distance_au for c in (_comet(s, giant_au=5.2) for s in range(400)) if c.period_class == "jupiter_family"]
    assert len(free_q) > 100 and len(held_q) > 100
    assert abs(sum(free_q) / len(free_q) - sum(held_q) / len(held_q)) > 0.2


def test_a_giant_no_orbit_of_the_family_can_be_held_by_gives_another_family():
    # A giant at 400 AU: T = 400 / a is far above 3 for every short orbit.
    for seed in range(60):
        comet = _comet(seed, giant_au=400.0)
        assert comet.period_class != "jupiter_family"
        assert comet.eccentricity >= 0.0


def test_the_ceiling_still_holds_under_a_giant():
    for seed in range(120):
        comet = _comet(seed, giant_au=5.2, ceiling_au=30.0)
        if comet.period_class == "jupiter_family":
            assert comet.perihelion_distance_au * (1 + comet.eccentricity) / (1 - comet.eccentricity) <= 30.0 + 1e-9


def test_the_redraw_kind_is_registered_and_the_tally_counts_it():
    assert "comet.tisserand" in KINDS
    with draw.counting() as tally:
        for seed in range(100):
            _comet(seed, giant_au=5.2)
    assert tally.get("comet.tisserand", 0) > 0
