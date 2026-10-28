"""
Generation knows the cloud a system sits in (GEN.357, `generation.cloud`): the field clouds that reach a sector are found
before it places a system, `SpaceSector.add_system` gives each system the cloud holding its position through
`StarSystem.apply_cloud`, and a system outside every cloud is changed in nothing.
"""

import json
import math

import pytest

from planetgen import tuning
from planetgen.db import query, store
from planetgen.galaxy import containment, nebula_field
from planetgen.galaxy.sector import SpaceSector
from planetgen.generation import cloud as cloudContext, fill_stages, run_galaxy, run_sector
from planetgen.generation.config import SystemConfig
from planetgen.generation.phenomena.nebula import Nebula
from planetgen.generation.system import StarSystem
from planetgen.physics import hii
from planetgen.physics.units import ly_to_pc, pc_to_ly
from planetgen.util import draw

from tests.test_db_persistence import _placed_nebula, _sector_with_one_system
from tests.test_nebula_field import ARM_ANGLE, EDGE_PC, REACH_PC, SEED, SHAPE, _on_ring

CENTER_PC = (0.0, 8000.0, 0.0)
"""tuple: A sector centre on the +y axis, so the sector's local +X is the galaxy's +y (`frame.sector_orientation`)."""


def _nebula(nebula_class, radius_ly=2.0):
    with draw.bound("gen357:nebula"):
        nebula = Nebula(SystemConfig(), nebula_class=nebula_class)
    nebula.radius_ly = radius_ly
    return nebula


def _system(planets=False):
    config = SystemConfig()
    config.PLANETS = planets
    with draw.bound("gen357:system"):
        return StarSystem(system_config=config)


def _sector(nebula, center_pc):
    """A filled-in sector (a fill has begun, so every stage has its stream) holding one field cloud at `center_pc`."""
    sector = SpaceSector("Gen357")
    with draw.bound("gen357:fill"):
        fill_stages.begin(sector)
    cloud_field = (SEED, SHAPE, CENTER_PC, REACH_PC)
    cloudContext.begin(sector, cloud_field, [(nebula, tuple(center_pc))])
    return sector


def _host_cloud(text, center=(0.0, 0.0, 0.0), cloud_class="D"):
    return cloudContext.cloud_dict("nebula", 7, "Test Nebula", cloud_class, "emission", 50.0, 9000.0, 0.5, 40.0, 3.0,
                                   center, text)


# --- The cloud dict -----------------------------------------------------------------------------------------------


def test_a_drawn_nebula_becomes_the_dict_a_loaded_system_holds():
    nebula = _nebula("D", radius_ly=30.0)
    cloud = cloudContext.cloud_of_nebula(nebula, (1.0, 2.0, 3.0))
    assert cloud["type"] == "nebula" and cloud["class"] == "D" and cloud["id"] is None and cloud["name"] is None
    assert cloud["radius_ly"] == 30.0 and cloud["age_myr"] == nebula.age_myr and cloud["center"] == (1.0, 2.0, 3.0)
    assert cloud["extinction_av"] == nebula.extinction_av and cloud["descriptor"] == nebula.nebula_type
    assert {star["type"] for star in cloud["host_stars"]} <= {name for name, _mass, _q in hii.TYPE_ANCHORS}
    assert sum(star["count"] for star in cloud["host_stars"]) > 0


def test_a_cloud_without_hosts_or_a_place_lists_no_host_stars():
    assert cloudContext.host_stars(None, (0.0, 0.0, 0.0)) == []
    assert cloudContext.host_stars("O5V x1", None) == []
    assert cloudContext.cloud_of_nebula(_nebula("M"), (0.0, 0.0, 0.0))["host_stars"] == []


def test_the_hosts_text_gives_each_type_its_count_luminosity_and_temperature():
    stars = cloudContext.host_stars("O5V x1, B0V x9, WN6 x1", (5.0, 6.0, 7.0))
    assert [(star["type"], star["count"]) for star in stars] == [("O5V", 1), ("B0V", 9)]   # a Wolf-Rayet type has no row
    assert stars[0]["luminosity_sol"] > stars[1]["luminosity_sol"] > 1.0e3
    assert stars[0]["teff_k"] > stars[1]["teff_k"] > 2.0e4
    assert all(star["position_pc"] == (5.0, 6.0, 7.0) for star in stars)


def test_a_remnant_row_gives_the_same_keys_as_a_nebula_row():
    remnant = cloudContext.cloud_dict("supernova_remnant", 3, "Old Shell", "S", "shell", 1.0, 2.0e6, 0.05, 20.0, 0.01,
                                      None)
    assert set(remnant) == set(cloudContext.cloud_of_nebula(_nebula("D"), (0, 0, 0)))
    assert remnant["center"] is None and remnant["host_stars"] == []


# --- G0, the tier and the enrichment ------------------------------------------------------------------------------


def test_the_far_ultraviolet_share_follows_the_design_note_and_is_held_beyond_its_ends():
    for teff, share in tuning.CLOUD_FUV_FRACTION_BY_TEFF_K.items():
        assert cloudContext.fuv_fraction(teff) == pytest.approx(share)
    assert cloudContext.fuv_fraction(5.0e3) == 0.16 and cloudContext.fuv_fraction(8.0e4) == 0.45
    assert cloudContext.fuv_fraction(25000.0) == pytest.approx((0.46 + 0.54) / 2)


def test_one_o9v_like_star_gives_about_a_thousand_habing_units_at_a_parsec():
    """The design note's table: 1e5 Lsun at 35,000 K is 1e3 G0 at 1 pc, falling as the square of the distance."""
    cloud = _host_cloud("O9V x1")
    cloud["host_stars"] = [{"type": "O9V", "count": 1, "luminosity_sol": 1.0e5, "teff_k": 35000.0,
                            "position_pc": (0.0, 0.0, 0.0)}]
    assert cloudContext.g0_at(cloud, (1.0, 0.0, 0.0)) == pytest.approx(1.0e3, rel=0.05)
    assert cloudContext.g0_at(cloud, (0.0, 10.0, 0.0)) == pytest.approx(10.0, rel=0.05)
    cloud["host_stars"][0]["count"] = 4
    assert cloudContext.g0_at(cloud, (1.0, 0.0, 0.0)) == pytest.approx(4.0e3, rel=0.05)


def test_a_host_is_never_taken_closer_than_the_minimum_distance():
    cloud = _host_cloud("O9V x1")
    at_centre = cloudContext.g0_at(cloud, (0.0, 0.0, 0.0))
    assert math.isfinite(at_centre) and at_centre == cloudContext.g0_at(
        cloud, (tuning.CLOUD_HOST_MIN_DISTANCE_PC / 2, 0.0, 0.0))


def test_the_g0_tier_counts_the_edges_reached():
    assert [cloudContext.g0_tier(g0) for g0 in (0.0, 99.0, 100.0, 999.0, 1.0e3, 1.0e4, 1.0e9)] == [0, 0, 1, 1, 2, 3, 3]


def test_the_tier_is_the_larger_of_the_class_tier_and_the_g0_tier():
    quiet = _host_cloud(None, cloud_class="R")
    assert cloudContext.effects_of(quiet, (0.0, 0.0, 0.0), draw).radiation_tier == tuning.CLOUD_CLASS_RADIATION_TIER["R"] == 3
    lit = _host_cloud("O9V x4, B0V x9", cloud_class="D")
    near = cloudContext.effects_of(lit, (0.5, 0.0, 0.0), draw)
    far = cloudContext.effects_of(lit, (500.0, 0.0, 0.0), draw)
    assert near.radiation_tier == cloudContext.g0_tier(near.g0) == 3
    assert far.g0 < 1.0 and far.radiation_tier == 0
    assert cloudContext.effects_of(_host_cloud(None, cloud_class="C"), (0.0, 0.0, 0.0), draw).radiation_tier == 1


def test_only_the_enriching_classes_draw_a_26al_factor():
    rng = draw.Stream("gen357:k")
    assert cloudContext.draw_enrichment("M", rng) is None and cloudContext.draw_enrichment("F", rng) is None
    assert cloudContext.draw_enrichment("D", rng) > 0.0


def test_the_enrichment_is_log_normal_about_the_field_median_with_the_classes_chance_of_the_rich_one():
    rng = draw.Stream("gen357:k2")
    plain = sorted(cloudContext.draw_enrichment("T", rng) for _ in range(4000))     # no rich chance for class T
    assert plain[2000] == pytest.approx(tuning.CLOUD_ENRICHMENT_FIELD_MEDIAN, rel=0.1)
    spread = [math.log10(k) for k in plain]
    mean = sum(spread) / len(spread)
    assert (sum((v - mean) ** 2 for v in spread) / len(spread)) ** 0.5 == pytest.approx(
        tuning.CLOUD_ENRICHMENT_SIGMA_DEX, rel=0.1)
    rich = sum(cloudContext.draw_enrichment("S", rng) > 0.3 for _ in range(4000)) / 4000   # chance 0.5 of the rich median
    assert 0.4 < rich < 0.6


# --- A system placed in a sector ----------------------------------------------------------------------------------


def test_a_system_placed_inside_a_field_cloud_gets_the_cloud_and_its_effects():
    nebula = _nebula("D")
    nebula.hosts = "O9V x2"
    centre = CENTER_PC
    sector = _sector(nebula, centre)
    system = _system()
    sector.add_system(system, position=(0.0, 0.0, 0.0))
    assert system.surrounding_cloud["class"] == "D" and system.surrounding_cloud["center"] == centre
    effects = system.cloud_effects
    assert effects.cloud_class == "D" and effects.cloud_id is None
    assert effects.enrichment_k is not None and effects.g0 > 0.0
    assert effects.radiation_tier == cloudContext.g0_tier(effects.g0)
    assert sector.stage_counts[(fill_stages.CLOUDS, "star-system")] == 1


def test_a_system_outside_every_cloud_is_untouched():
    sector = _sector(_nebula("D"), CENTER_PC)
    system = _system()
    sector.add_system(system, position=(8.0, 0.0, 0.0))        # 8 ly out, against a 2 ly cloud
    assert system.surrounding_cloud is None and system.cloud_effects is None
    assert (fill_stages.CLOUDS, "star-system") not in sector.stage_counts


def test_the_sectors_local_axes_are_turned_into_the_galaxy_frame_before_the_test():
    """The cloud lies 3 ly along the sector's local +X, which is the galaxy's +y for a centre on the +y axis."""
    nebula = _nebula("D", radius_ly=2.0)
    sector = _sector(nebula, (0.0, CENTER_PC[1] + ly_to_pc(3.0), 0.0))
    along_x, along_y = _system(), _system()
    sector.add_system(along_x, position=(3.0, 0.0, 0.0))
    sector.add_system(along_y, position=(0.0, 3.0, 0.0))
    assert along_x.surrounding_cloud is not None
    assert along_y.surrounding_cloud is None


def test_the_smaller_of_two_nested_clouds_is_the_one_a_system_gets():
    big, small = _nebula("D", radius_ly=30.0), _nebula("C", radius_ly=2.0)
    sector = SpaceSector("Nested")
    with draw.bound("gen357:fill"):
        fill_stages.begin(sector)
    cloudContext.begin(sector, (SEED, SHAPE, CENTER_PC, REACH_PC), [(big, CENTER_PC), (small, CENTER_PC)])
    system = _system()
    sector.add_system(system, position=(0.0, 0.0, 0.0))
    assert system.surrounding_cloud["class"] == "C"


def test_a_cloud_draws_the_enrichment_from_a_stream_of_its_own():
    """Placing a system inside a cloud moves no draw of the stream bound at the time."""
    sector = _sector(_nebula("D"), CENTER_PC)
    system = _system()
    with draw.bound("gen357:bound"):
        sector.add_system(system, position=(0.0, 0.0, 0.0))
        after = [draw.random() for _ in range(3)]
    with draw.bound("gen357:bound"):
        assert [draw.random() for _ in range(3)] == after


def test_the_same_placement_draws_the_same_enrichment_and_another_placement_another():
    def enrichment(fill_seed, entries=0):
        sector = SpaceSector("Repeat")
        with draw.bound(fill_seed):
            fill_stages.begin(sector)
        cloudContext.begin(sector, (SEED, SHAPE, CENTER_PC, REACH_PC), [(_nebula("D"), CENTER_PC)])
        system = _system()
        sector.add_system(system, position=(0.0, 0.0, 0.0))
        return system.cloud_effects.enrichment_k

    assert enrichment("a") == enrichment("a")
    assert enrichment("a") != enrichment("b")


def test_a_sector_without_a_cloud_field_has_no_containers_and_no_work_in_add_system():
    sector = SpaceSector("Plain")
    assert sector.cloud_containers == []
    system = _system()
    sector.add_system(system, position=(0.0, 0.0, 0.0))
    assert system.surrounding_cloud is None and system.cloud_effects is None


def test_apply_cloud_leaves_the_planets_as_drawn():
    system = _system(planets=True)
    before = json.dumps(system.to_dict(), sort_keys=True, default=str)
    system.apply_cloud(cloudContext.cloud_of_nebula(_nebula("D"), CENTER_PC), CENTER_PC, rng=draw.Stream("gen357:keep"))
    assert json.dumps(system.to_dict(), sort_keys=True, default=str) == before


# --- A whole sector ------------------------------------------------------------------------------------------------


def _generate(monkeypatch, cloud_in_play, hook=True):
    """A galaxy-placed sector of planets-and-all systems with one field cloud over its centre (or none), the hook on or off."""
    from planetgen.galaxy.spherical import SectorCell, sector_address_at, sector_position_pc

    address = sector_address_at(_on_ring(ARM_ANGLE), EDGE_PC)
    position = sector_position_pc(*address, EDGE_PC)
    clouds = [(_nebula("D", radius_ly=76.0), tuple(position))] if cloud_in_play else []
    monkeypatch.setattr(nebula_field, "clouds_reaching", lambda *args: clouds)
    if not hook:
        monkeypatch.setattr(cloudContext, "apply_in_sector", lambda *args: None)
    args = run_galaxy._default_generation_args()
    args.num_systems = 6
    args.density = None
    cell = SectorCell.for_address(*address, pc_to_ly(EDGE_PC))
    with draw.bound("gen357:sector"):
        return run_sector.generate_sector(args, 25000.0, cell=cell,
                                          cloud_field=(SEED, SHAPE, tuple(position), REACH_PC))[1]


def test_a_generated_sector_inside_a_seeded_cloud_hands_every_system_the_cloud(monkeypatch):
    sector = _generate(monkeypatch, True)
    assert sector.entries
    assert all(entry.star_system.surrounding_cloud["class"] == "D" for entry in sector.entries)
    assert sector.stage_counts[(fill_stages.CLOUDS, "star-system")] == len(sector.entries)
    assert all(entry.star_system.cloud_effects.enrichment_k is not None for entry in sector.entries)


def test_the_cloud_changes_nothing_a_system_is_made_of(monkeypatch):
    with_cloud = _generate(monkeypatch, True)
    without_hook = _generate(monkeypatch, True, hook=False)
    assert all(entry.star_system.surrounding_cloud is None for entry in without_hook.entries)
    assert (json.dumps(with_cloud.to_dict(), sort_keys=True, default=str)
            == json.dumps(without_hook.to_dict(), sort_keys=True, default=str))


def test_a_sector_no_cloud_reaches_has_no_containers_and_no_cloud_stage(monkeypatch):
    sector = _generate(monkeypatch, False)
    assert sector.cloud_containers == [] and sector.field_nebulae == []
    assert all(entry.star_system.surrounding_cloud is None for entry in sector.entries)
    assert not any(stage == fill_stages.CLOUDS for stage, _kind in sector.stage_counts)


def test_the_clouds_are_found_once_per_sector(monkeypatch):
    calls = []
    real = nebula_field.clouds_reaching
    monkeypatch.setattr(nebula_field, "clouds_reaching", lambda *args: calls.append(args) or real(*args))
    from planetgen.galaxy.spherical import SectorCell, sector_address_at, sector_position_pc
    address = sector_address_at(_on_ring(ARM_ANGLE), EDGE_PC)
    position = sector_position_pc(*address, EDGE_PC)
    args = run_galaxy._default_generation_args()
    args.num_systems, args.density = 2, None
    with draw.bound("gen357:once"):
        run_sector.generate_sector(args, 25000.0, cell=SectorCell.for_address(*address, pc_to_ly(EDGE_PC)),
                                   cloud_field=(SEED, SHAPE, tuple(position), REACH_PC))
    assert len(calls) == 1


# --- The store -----------------------------------------------------------------------------------------------------


def test_a_loaded_system_holds_the_extended_cloud_and_the_listing_keeps_its_old_shape(mysql_config):
    sector_id = store.save_sector(_sector_with_one_system("Inside"), config=mysql_config, galaxy_position={
        "center_x_pc": 100.0, "center_y_pc": 0.0, "center_z_pc": 0.0, "galactic_radius_pc": 100.0,
    })
    conn = store.get_connection(mysql_config)
    try:
        system_id = conn.execute("SELECT id FROM star_systems WHERE sector_id = ?", (sector_id,)).fetchone()["id"]
        with conn:
            nebula_id = _placed_nebula(conn, sector_id, (100.0, 0.0, 0.0), radius_ly=50.0)
            conn.execute("UPDATE nebulae SET extinction_av = 2.5, age_myr = 4.0, hosts = 'O9V x2' WHERE id = ?",
                         (nebula_id,))
        cloud = store.load_star_system(conn, system_id).surrounding_cloud
        assert cloud["id"] == nebula_id and cloud["extinction_av"] == 2.5 and cloud["age_myr"] == 4.0
        row = conn.execute("SELECT center_x_pc, center_y_pc, center_z_pc FROM nebulae WHERE id = ?", (nebula_id,)).fetchone()
        assert cloud["radius_ly"] == 50.0
        assert cloud["center"] == (row["center_x_pc"], row["center_y_pc"], row["center_z_pc"])
        assert [(star["type"], star["count"]) for star in cloud["host_stars"]] == [("O9V", 2)]
        listed = query.containing_cloud(conn, {"inside_nebula_id": nebula_id, "inside_remnant_id": None})
        assert tuple(listed) == cloudContext.LISTED_KEYS
    finally:
        conn.close()


def test_containment_is_one_function_for_the_store_and_for_generation():
    near = {"column": "inside_nebula_id", "id": 1, "center": (0.0, 0.0, 0.0), "radius_pc": 4.0}
    assert containment.innermost_container((1.0, 0.0, 0.0), [near]) is near
    assert not hasattr(store, "innermost_container")
