# tests/test_belt_roster.py

"""
GEN.1070: where each stored member of a belt is (`physics/belt_roster.py`) and the members a belt keeps (`AsteroidBelt.members`).
"""

import math

import pytest

from planetgen.generation.belt import AsteroidBelt
from planetgen.physics import belt_ledger, belt_roster, constants
from planetgen.util import draw

from tests.fuzz_support import deterministic_entropy
from tests.test_belt_birth import _belts, _config
from planetgen.generation.system import StarSystem

STAR = constants.SOLAR_MASS_TO_KG
SNOW = 2.7


def _roster(count=50, lower=2.0, upper=3.4, key="a"):
    return belt_roster.draw_roster(count, lower, upper, STAR, SNOW, 1.0, draw.Stream(key))


def test_every_body_lies_inside_the_belt_with_bounded_elements():
    for entry in _roster(count=300, lower=2.0, upper=6.0):
        assert 2.0 <= entry.placement.distance_au <= 6.0
        assert 0.0 <= entry.eccentricity <= belt_roster.E_MAX
        assert 0.0 <= entry.inclination_deg <= belt_roster.I_MAX_DEG
        assert all(0.0 <= angle < 360.0 for angle in (entry.node_deg, entry.periapsis_deg, entry.anomaly_deg))
        assert 0.0 <= entry.placement.u_taxon < 1.0 and 0.0 <= entry.placement.u_porosity < 1.0


def test_the_same_belt_gives_the_same_roster_and_a_shorter_ledger_the_same_first_bodies():
    assert _roster(300) == _roster(300)
    assert _roster(300)[:20] == _roster(20)
    assert _roster(10, key="a") != _roster(10, key="b")


def test_the_draws_taken_do_not_depend_on_the_count():
    one, two = draw.Stream("same"), draw.Stream("same")
    belt_roster.draw_roster(3, 2.0, 3.0, STAR, SNOW, 1.0, one)
    belt_roster.draw_roster(300, 2.0, 3.0, STAR, SNOW, 1.0, two)
    assert one.random() == two.random()


def test_a_wider_belt_has_a_larger_eccentricity_and_inclination_on_the_average():
    def mean(entries, name):
        return math.fsum(getattr(entry, name) for entry in entries) / len(entries)
    narrow, wide = _roster(300, 2.5, 2.7), _roster(300, 2.0, 6.0)
    assert mean(wide, "eccentricity") > 3.0 * mean(narrow, "eccentricity")
    assert mean(wide, "inclination_deg") > 3.0 * mean(narrow, "inclination_deg")


def test_the_distances_follow_where_the_solids_are():
    entries = _roster(300, 2.0, 6.0)
    inner = sum(entry.placement.distance_au < 3.0 for entry in entries)
    table = belt_roster.cumulative_table(2.0, 6.0, STAR, SNOW, 1.0)
    share = belt_roster.semi_major_axis_au(table, 0.5)
    assert 2.0 < share < 6.0
    expected = sum(1 for u in (i / 300 for i in range(300)) if belt_roster.semi_major_axis_au(table, u) < 3.0)
    assert abs(inner - expected) < 60       # a draw of 300 numbers from the same distribution


def test_the_table_is_inverted_at_its_ends():
    table = belt_roster.cumulative_table(2.0, 3.4, STAR, SNOW, 1.0)
    assert belt_roster.semi_major_axis_au(table, 0.0) == pytest.approx(2.0)
    assert belt_roster.semi_major_axis_au(table, 1.0 - 1e-12) == pytest.approx(3.4, rel=1e-6)


def test_a_bad_belt_or_count_is_refused():
    with pytest.raises(belt_roster.BeltRosterError):
        belt_roster.draw_roster(5, 3.0, 2.0, STAR, SNOW, 1.0, draw.Stream("x"))
    with pytest.raises(belt_roster.BeltRosterError):
        belt_roster.draw_roster(belt_ledger.LEDGER_MAX_BODIES + 1, 2.0, 3.0, STAR, SNOW, 1.0, draw.Stream("x"))


def _generated_belt():
    config = _config()
    with deterministic_entropy(1070):
        for seed in range(1, 40):
            draw.set_run_seed(seed)
            system = StarSystem(system_config=config)
            belt = next((b for b in _belts(system) if b.origin is not None), None)
            if belt is not None:
                return belt
    pytest.fail("no belt with an origin in 39 systems")


def test_a_generated_belt_keeps_members_heaviest_first_inside_its_limits():
    belt = _generated_belt()
    members = belt.members
    assert members and len(members) <= belt_ledger.LEDGER_MAX_BODIES
    ranks = [member["rank"] for member in members]          # GEN.1510, GEN.1693: a body that became a planet keeps its rank and leaves a gap (one inside its ring and clear of the other plane is not always the heaviest), so the ranks left rise without a repeat
    assert ranks[0] >= 1 and ranks == sorted(set(ranks))
    masses = [member["mass_kg"] for member in members]
    assert masses == sorted(masses, reverse=True)
    assert sum(masses) < belt.mass_kg * 1.0000001
    for member in members:
        assert belt.lower_limit <= member["distance_au"] <= belt.upper_limit
        assert member["radius_km"] > 0.0 and 0.0 < member["density_g_cm3"] < 10.0
        shares = member["rock_fraction"] + member["metal_fraction"] + member["water_fraction"] + member["ice_fraction"]
        assert shares == pytest.approx(1.0, abs=1e-6)
        assert member["status"] in ("planet", "dwarf planet", "minor body")


def test_a_belt_cut_back_draws_its_members_again_for_the_shorter_belt():
    belt = _generated_belt()
    before = belt.members
    ceiling = (belt.lower_limit * belt.upper_limit) ** 0.5
    if belt.shorten_to(ceiling):
        after = belt.members
        assert all(member["distance_au"] <= ceiling for member in after)
        assert after is not before


def test_a_belt_read_back_keeps_the_members_it_was_saved_with():
    belt = _generated_belt()
    members = belt.members
    again = AsteroidBelt.from_dict(belt.to_dict(), belt.system_config)
    assert again.members == members
    bare = belt.to_dict()
    bare["members"] = []
    assert AsteroidBelt.from_dict(bare, belt.system_config).members == []


# --- GEN.1092: the members are not placed in what the giants cleared -------------------------------------------------------------------------

def test_no_body_is_placed_in_a_cleared_stretch_and_the_rest_keep_their_share_of_the_solids():
    table = belt_roster.cumulative_table(2.0, 4.0, 1.989e30, 2.7, 1.0)
    cleared = [(2.4, 2.6), (3.1, 3.3)]
    places = [belt_roster.semi_major_axis_au(table, (k + 0.5) / 2000.0, cleared) for k in range(2000)]
    assert not any(low < a < high for a in places for low, high in cleared)
    assert places == sorted(places) and 2.0 <= places[0] and places[-1] <= 4.0
    free = belt_roster.mass_below(table, 4.0) - sum(belt_roster.mass_below(table, high) - belt_roster.mass_below(table, low) for low, high in cleared)
    inner = sum(1 for a in places if a < 2.4) / 2000.0
    assert inner == pytest.approx(belt_roster.mass_below(table, 2.4) / free, abs=0.01)


def test_with_nothing_cleared_the_table_is_read_as_before():
    table = belt_roster.cumulative_table(2.0, 4.0, 1.989e30, 2.7, 1.0)
    assert belt_roster.semi_major_axis_au(table, 0.37, ()) == belt_roster.semi_major_axis_au(table, 0.37)
    assert belt_roster.mass_below(table, 2.0) == 0.0 and belt_roster.mass_below(table, 9.0) == table[1][-1]
    for u in (0.1, 0.5, 0.9):
        assert belt_roster.mass_below(table, belt_roster.semi_major_axis_au(table, u)) == pytest.approx(u * table[1][-1], rel=1e-9)
