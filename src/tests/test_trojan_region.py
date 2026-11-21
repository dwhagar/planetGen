# tests/test_trojan_region.py

"""
GEN.1150: the trojan swarm of a giant is made into the regions at its L4 and L5 (`generation/trojan_region.py`), the mass from the outer disk the giant swept (`swept_hand_on.trojan_part`), stored as belt rows
whose parent is the giant, and named for it (`<host> L4 trojan region`).
"""

import math
import types

import pytest

from planetgen.db import store
from planetgen.generation import belt_birth, trojan_region
from planetgen.generation.belt import AsteroidBelt
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.generation.trojan_region import TrojanRegion
from planetgen.names import bodies as body_names
from planetgen.physics import scatter, swept_hand_on, trojan_swarm, trojans
from planetgen.physics import unit_registry as units
from planetgen.util import draw
from tests.test_systems import make_config

JUPITER_KG = 1.898e27


@pytest.fixture(scope="module")
def system():
    with draw.bound(5):
        return StarSystem(system_config=make_config("G2V", PLANETS=False))


def _giant(au=5.2, mass_kg=JUPITER_KG):
    return types.SimpleNamespace(body_type="g", distance=au, mass=mass_kg, name="Test III")


def _disk(trojan_kg, giant):
    return types.SimpleNamespace(body_type="a", swarm_kg=trojan_kg, scatter_giant_kg=giant.mass, scatter_giant_au=giant.distance)


def test_the_swarm_is_the_bound_share_of_what_the_giant_swept_times_the_trojan_row_and_the_neighbours_share():
    theta = 20.0
    ejected = swept_hand_on.ejected_share(theta)
    whole = swept_hand_on.trojan_part(1.0e24, theta)
    assert whole == pytest.approx(1.0e24 * (1.0 - ejected) * dict(scatter.SWEPT_SPLIT)["trojans"])
    assert swept_hand_on.trojan_part(1.0e24, theta, 0.25) == pytest.approx(0.25 * whole)
    assert swept_hand_on.trojan_part(1.0e24, theta, 7.0) == pytest.approx(whole) and swept_hand_on.trojan_part(1.0e24, theta, -1.0) == 0.0     # the share is kept between 0 and 1
    assert swept_hand_on.trojan_part(0.0, theta) == 0.0 and swept_hand_on.trojan_part(1.0e24, None) == 0.0 and swept_hand_on.trojan_part(1.0e24, 0.0) == 0.0


def test_the_swarm_is_inside_other_and_the_scatter_step_gives_the_same_swarm_from_the_swept_share():
    theta = 20.0
    parts = swept_hand_on.hand_on(1.0e24, theta)
    assert swept_hand_on.trojan_part(1.0e24, theta) < parts["other"]
    shares = scatter.scatter_shares(1.0, 300.0, 5.2, 11.0, 0.0, 1.0, 100.0)
    swept = shares["ejected"] + sum(shares[key] for key, _ in scatter.SWEPT_SPLIT)
    ejected = swept_hand_on.ejected_share(scatter.safronov_number(300.0, 11.0, 5.2, 1.0))
    assert shares["trojans"] / swept == pytest.approx((1.0 - ejected) * dict(scatter.SWEPT_SPLIT)["trojans"])


def test_a_giant_with_a_swept_swarm_gets_a_leading_and_a_trailing_region_about_its_orbit(system):
    giant = _giant()
    found = trojan_region.build([giant, _disk(1.0e22, giant)], system.star, system.system_config)
    assert [region.lagrange_point for region in found] == ["L4", "L5"]
    for region in found:
        assert isinstance(region, TrojanRegion) and region.body_type == "a" and region.host is giant
        assert region.distance == giant.distance
        assert region.lower_limit == pytest.approx(giant.distance * (1.0 - region.mean_eccentricity)) and region.upper_limit == pytest.approx(giant.distance * (1.0 + region.mean_eccentricity))
        assert 0.0 < region.mass_kg <= 1.0e22 and region.largest_body_km > 0.0 and region.composition and region.origin == trojan_region.ORIGIN
        assert 0.0 <= region.libration_amplitude_deg <= trojans.TADPOLE_LIMIT_DEG and 0.0 <= region.libration_phase_deg < 360.0
        assert 0.0 <= region.mean_eccentricity <= trojans.ECC_MAX and 0.0 <= region.mean_inclination_deg <= trojans.INCLINATION_MAX_DEG
    assert sum(region.primordial_mass_earth for region in found) == pytest.approx(units.convert(1.0e22, "kg", "Mearth"))


def test_the_regions_are_drawn_from_the_swarms_own_stream_so_the_same_giant_gets_the_same_ones(system):
    giant = _giant()
    first = trojan_region.build([giant, _disk(1.0e22, giant)], system.star, system.system_config)
    draw.random()                  # a draw elsewhere moves nothing of the swarm
    second = trojan_region.build([giant, _disk(1.0e22, giant)], system.star, system.system_config)
    assert [r.to_dict() for r in first] == [r.to_dict() for r in second]


def test_no_swarm_a_swarm_too_light_to_keep_a_pair_that_cannot_hold_it_and_a_giant_that_swept_nothing_give_no_region(system):
    giant = _giant()
    assert trojan_region.build([giant, _disk(0.0, giant)], system.star, system.system_config) == []
    assert trojan_region.build([giant], system.star, system.system_config) == []
    assert trojan_region.build([giant, _disk(1.0e5, giant)], system.star, system.system_config) == []                    # under the belt minimum once ground
    heavy = _giant(mass_kg=0.05 * system.star.mass)                      # a secondary over 4 percent of the star holds nothing at L4 and L5
    assert trojan_region.build([heavy, _disk(1.0e22, heavy)], system.star, system.system_config) == []


def test_the_outermost_giant_holds_the_swarm(system):
    inner, outer = _giant(2.0), _giant(9.0)
    found = trojan_region.build([inner, outer, _disk(1.0e22, outer)], system.star, system.system_config)
    assert found and all(region.host is outer for region in found)


def test_a_region_read_back_from_its_dict_keeps_everything_it_had(system):
    giant = _giant()
    region = trojan_region.build([giant, _disk(1.0e22, giant)], system.star, system.system_config)[0]
    again = TrojanRegion.from_dict(region.to_dict(), system.system_config)
    assert again.to_dict() == region.to_dict() and again.members == []


def test_a_region_has_no_members_and_a_plain_belt_has_none_of_its_fields(system):
    giant = _giant()
    region = trojan_region.build([giant, _disk(1.0e22, giant)], system.star, system.system_config)[0]
    assert region.members == []
    with draw.bound(5):
        belt = AsteroidBelt(system.system_config, system.star, 2.5, 2.2, 3.4)
    assert belt.lagrange_point is None and belt.libration_amplitude_deg is None and belt.mean_eccentricity is None and "lagrange_point" not in belt.to_dict()


def test_a_region_is_named_for_its_giant_and_takes_no_number_among_the_belts():
    assert body_names.region_label("Voranthis V", "L4") == "Voranthis V L4 trojan region"
    belts = [{"id": 1, "star_id": 7, "orbital_index": 2, "lagrange_point": None}, {"id": 2, "star_id": 7, "orbital_index": 3, "lagrange_point": "L4"},
             {"id": 3, "star_id": 7, "orbital_index": 3, "lagrange_point": "L5"}, {"id": 4, "star_id": 7, "orbital_index": 6, "lagrange_point": None}]
    assert [body_names.belt_number(belt, belts) for belt in (belts[0], belts[3])] == [1, 2]


def test_the_survival_of_the_swarm_is_the_one_the_scatter_step_uses():
    giants = [types.SimpleNamespace(mass=JUPITER_KG, distance=5.2), types.SimpleNamespace(mass=JUPITER_KG, distance=5.6)]
    star = types.SimpleNamespace(mass=1.989e30)
    assert belt_birth._swarm_survival(star, star.mass, giants) == trojans.surviving_share(star.mass, JUPITER_KG, 5.6, [(JUPITER_KG, 5.2)], None) < 1.0
    assert belt_birth._swarm_survival(star, star.mass, giants[:1]) == 1.0


def test_a_region_is_no_belt_of_the_count_a_stored_moon_or_a_promoted_planet_names_a_belt_by():
    from planetgen.generation import irregular_moons
    plain = [types.SimpleNamespace(body_type="a", lagrange_point=None), types.SimpleNamespace(body_type="a")]
    region = types.SimpleNamespace(body_type="a", lagrange_point="L4")
    found = irregular_moons.belts_of([plain[0], region, plain[1]])
    assert list(found) == [0, 1] and found[0] is plain[0] and found[1] is plain[1]


@pytest.mark.fresh_database
def test_regions_are_saved_with_their_giant_as_parent_and_read_back_beside_it(mysql_config):
    config = SystemConfig()
    config.PLANETS = True
    config.MAX_PLANETS = True
    config.BINARY_SYSTEM = False
    config.FORCE_CLASS = "J"
    with draw.bound(3):
        made = StarSystem(system_config=config)
    giant = next(planet for planet in made.primary.orbiters if planet.body_type == "g")
    regions = [TrojanRegion(config, made.star, giant, region)
               for region in trojan_swarm.draw_regions(made.star.mass, giant.mass, 1.0e22, draw.Stream("test:trojan-store"))]
    assert len(regions) == 2
    made.primary.orbiters.extend(regions)
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            system_id = store.insert_star_system(conn, made, config)
        rows = conn.execute("SELECT * FROM asteroid_belts WHERE star_system_id = ? AND lagrange_point IS NOT NULL ORDER BY lagrange_point", (system_id,)).fetchall()
        giant_row = conn.execute("SELECT id FROM planets WHERE star_system_id = ? AND orbital_index = ?", (system_id, made.primary.orbiters.index(giant))).fetchone()
        assert [row["lagrange_point"] for row in rows] == ["L4", "L5"]
        assert all(row["parent_type"] == "planet" and row["parent_id"] == giant_row["id"] for row in rows)
        assert [row["libration_amplitude_deg"] for row in rows] == [region.libration_amplitude_deg for region in regions]
        loaded = store.load_star_system(conn, system_id)
        back = [body for body in loaded.primary.orbiters if getattr(body, "lagrange_point", None)]
        assert [body.lagrange_point for body in back] == ["L4", "L5"] and all(isinstance(body, TrojanRegion) for body in back)
        assert all(body.host is not None and body.host.body_type == "g" and body.host.db_id == giant_row["id"] for body in back)
        assert [body.mass_kg for body in back] == [region.mass_kg for region in regions] and [body.mean_eccentricity for body in back] == [region.mean_eccentricity for region in regions]
        plain = [body for body in loaded.primary.orbiters if body.body_type == "a" and not getattr(body, "lagrange_point", None)]
        assert all(type(body) is AsteroidBelt for body in plain)
        from planetgen.db import belt_names
        names = belt_names.names_for_belts(conn, [row["id"] for row in rows])
        assert sorted(names.values()) == [f"{giant.name} L4 trojan region", f"{giant.name} L5 trojan region"]
    finally:
        conn.rollback()
        conn.close()


def test_a_trojan_region_and_a_trojan_planet_share_their_giants_orbit_and_are_no_slot_of_the_spacing_rules():
    """GEN.1150, GEN.1151: a region is put at the end of its star's list and stands at its giant's distance, and a trojan planet stands at it too: `check_orbit_spacing` must not call either out of order or too close."""
    from planetgen.generation import validation
    outer = types.SimpleNamespace(body_type="a", distance=24.13, lower_limit=24.13, upper_limit=30.0, mass=0.0, lagrange_point=None, belt_rank=None, radius=0)
    region = types.SimpleNamespace(body_type="a", distance=15.0, lower_limit=14.85, upper_limit=15.15, mass=0.0, lagrange_point="L4", belt_rank=None, radius=0)
    assert validation.check_orbit_spacing([outer, region]) == []
    assert validation.slot_planets([outer, region]) == [outer]
    region.lagrange_point = None         # the same body with no point is a belt listed out of order: the check still sees it
    assert [problem.message for problem in validation.check_orbit_spacing([outer, region])] == ["it is listed after asteroid belt at 24.13 AU but orbits inside it"]
