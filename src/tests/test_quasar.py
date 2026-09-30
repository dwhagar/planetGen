# tests/test_quasar.py

"""
Tests for the quasar phenomenon (schema v31): the `Quasar` model itself,
its one-per-galaxy placement at the galactic center, storage, and how it
shows up in the sector's phenomena and on the site.
"""

import math

import pytest

import generate
import queryDb
from stellarObjects import _db, program_constants
from stellarObjects.config import SystemConfig
from stellarObjects.galaxyGeometry import sector_position_pc, shell_sector_count
from stellarObjects.quasarData import Quasar
from stellarObjects.spaceSector import SpaceSector
from stellarObjects.utils import pc_to_ly
from tests.test_galaxy_gen import EDGE_PC, _seed_skeleton


def _core_sector(mysql_config, slot=0, sector=None):
    """Saves `sector` (empty by default) at a real shell-0 address, the
    way `generate_and_save_sector_at` would."""
    x, y, z = sector_position_pc(0, slot, EDGE_PC)
    galaxy_position = {
        "center_x_pc": x, "center_y_pc": y, "center_z_pc": z,
        "galactic_radius_pc": math.sqrt(x * x + y * y + z * z),
        "shell_index": 0, "shell_slot_index": slot,
        "vertices_pc": {"inner": [], "outer": []},
    }
    if sector is None:
        sector = SpaceSector(f"Core {slot}")
    return _db.save_sector(sector, config=mysql_config, galaxy_position=galaxy_position)


def test_quasar_physics_follows_from_mass_and_eddington_ratio():
    for _ in range(50):
        q = Quasar(SystemConfig())
        lo, hi = program_constants.QUASAR_BLACK_HOLE_MASS_RANGE_SOLAR
        assert lo <= q.black_hole_mass_solar <= hi
        assert q.luminosity_w == pytest.approx(
            q.eddington_ratio * program_constants.EDDINGTON_LUMINOSITY_W_PER_SOLAR_MASS * q.black_hole_mass_solar
        )
        # Every quasar outshines its host galaxy's starlight.
        assert q.galaxy_luminosity_multiple > 1
        assert (q.jet_length_ly is not None) == q.is_radio_loud
        assert q.broad_line_region_light_days > 0


def test_quasar_round_trips_and_describes_itself():
    cfg = SystemConfig()
    q = Quasar(cfg, name="Test Nucleus")
    copy = Quasar.from_dict(q.to_dict(), cfg)
    assert copy.to_dict() == q.to_dict()
    text = str(q)
    assert "Test Nucleus (Quasar)" in text
    assert "active nucleus" in text


def test_space_sector_serialization_keeps_a_quasar():
    sector = SpaceSector("Round Trip")
    sector.add_phenomenon(Quasar(SystemConfig()), "quasar", position=(0.0, 0.0, -1.0))
    rebuilt = SpaceSector.from_dict(sector.to_dict())
    assert isinstance(rebuilt.phenomena[0].phenomenon, Quasar)


def test_random_phenomenon_choice_never_picks_a_quasar():
    assert "quasar" in program_constants.PHENOMENON_TYPE_CHOICES
    assert "quasar" not in program_constants.RANDOM_PHENOMENON_TYPE_CHOICES
    assert "quasar" not in program_constants.PHENOMENON_RATE_PER_STAR_SYSTEM


def test_add_galactic_nucleus_respects_the_chance(monkeypatch):
    args = generate._default_generation_args()
    monkeypatch.setattr(program_constants, "QUASAR_ACTIVE_NUCLEUS_CHANCE", 0.0)
    sector = SpaceSector("Quiet Core")
    assert generate.add_galactic_nucleus(sector, args, 4.0) is None
    assert sector.phenomena == []

    monkeypatch.setattr(program_constants, "QUASAR_ACTIVE_NUCLEUS_CHANCE", 1.0)
    entry = generate.add_galactic_nucleus(sector, args, 4.0)
    assert entry.phenomenon_type == "quasar"
    assert entry.position == (0.0, 0.0, -4.0)


def test_core_sector_quasar_is_stored_at_the_galactic_center(mysql_config, monkeypatch):
    # The in-sector offset add_galactic_nucleus picks must convert back to
    # the galactic origin through the sector's own rotated cube frame.
    x, y, z = sector_position_pc(0, 0, EDGE_PC)
    distance_ly = pc_to_ly(math.sqrt(x * x + y * y + z * z))

    monkeypatch.setattr(program_constants, "QUASAR_ACTIVE_NUCLEUS_CHANCE", 1.0)
    sector = SpaceSector("Core Host")
    generate.add_galactic_nucleus(sector, generate._default_generation_args(), distance_ly)
    placement = _db._galaxy_placement_from_sector_offset(
        {"center_x_pc": x, "center_y_pc": y, "center_z_pc": z}, sector.phenomena[0].position,
    )
    assert placement["galactic_radius_pc"] == pytest.approx(0.0, abs=1e-9)

    sector_id = _core_sector(mysql_config, sector=sector)
    conn = _db.get_connection(mysql_config)
    try:
        row = conn.execute("SELECT * FROM quasars WHERE sector_id = ?", (sector_id,)).fetchone()
        matches = queryDb.phenomena_near_sector(conn, sector_id)
        detail = queryDb.phenomenon_detail(conn, "quasar", row["id"])
    finally:
        conn.close()

    assert (row["center_x_pc"], row["center_y_pc"], row["center_z_pc"], row["galactic_radius_pc"]) == (0, 0, 0, 0)
    assert [(m["type"], m["id"]) for m in matches] == [("quasar", row["id"])]
    assert matches[0]["descriptor"] in ("radio-loud", "radio-quiet")
    assert detail["black_hole_mass_solar"] == row["black_hole_mass_solar"]


def test_only_the_first_core_sector_rolls_for_a_quasar(mysql_config, monkeypatch):
    n_0 = shell_sector_count(0)
    _seed_skeleton(mysql_config, bands=[(0, 0, 0, n_0 - 1)])
    monkeypatch.setattr(program_constants, "QUASAR_ACTIVE_NUCLEUS_CHANCE", 1.0)
    monkeypatch.setattr(
        generate, "generate_sector",
        lambda args, galactic_center_dist_ly=None: ("Fake", SpaceSector("Fake")),
    )

    for slot in range(n_0):
        generate.ensure_sector_generated(0, slot, config=mysql_config)

    conn = _db.get_connection(mysql_config)
    try:
        rows = conn.execute(
            "SELECT q.galactic_radius_pc, s.shell_slot_index FROM quasars q JOIN sectors s ON s.id = q.sector_id"
        ).fetchall()
    finally:
        conn.close()
    assert len(rows) == 1
    assert rows[0]["shell_slot_index"] == 0
    assert rows[0]["galactic_radius_pc"] == pytest.approx(0.0, abs=1e-9)


def test_save_phenomenon_places_a_quasar_only_at_the_core_and_only_once(mysql_config):
    cfg = SystemConfig()
    outer_id = _db.save_sector(
        SpaceSector("Outer"), config=mysql_config,
        galaxy_position={
            "center_x_pc": 50.0, "center_y_pc": 0.0, "center_z_pc": 0.0, "galactic_radius_pc": 50.0,
            "shell_index": 5, "shell_slot_index": 0, "vertices_pc": {"inner": [], "outer": []},
        },
    )
    with pytest.raises(ValueError, match="shell-0"):
        _db.save_phenomenon(Quasar(cfg), cfg, "quasar", config=mysql_config, sector_id=outer_id)

    core_id = _core_sector(mysql_config, slot=1)
    quasar_id = _db.save_phenomenon(Quasar(cfg), cfg, "quasar", config=mysql_config, sector_id=core_id)
    with pytest.raises(ValueError, match="already has a quasar"):
        _db.save_phenomenon(Quasar(cfg), cfg, "quasar", config=mysql_config, sector_id=core_id)

    # Unplaced quasars (no sector) are still allowed, like every other type.
    _db.save_phenomenon(Quasar(cfg), cfg, "quasar", config=mysql_config)

    conn = _db.get_connection(mysql_config)
    try:
        row = conn.execute("SELECT * FROM quasars WHERE id = ?", (quasar_id,)).fetchone()
        assert queryDb.count_phenomena(conn) == 2
    finally:
        conn.close()
    # Snapped to the center rather than jittered around the sector.
    assert (row["center_x_pc"], row["center_y_pc"], row["center_z_pc"]) == (0, 0, 0)

