# tests/test_hosted_nebula_rows.py

"""
GEN.530: the nebula a pre-placed star lights shows from its `bright_stars` row. The row rolls the nebula its sector
builds (`nebula_contents.row_nebula`), the star's page and the uncharted sector's scene read it before any fill, and a
star already built into a system lists none (the system's own phenomenon holds it).
"""

import types

import numpy as np
import pytest

from planetgen.db import hostednebula, query, store
from planetgen.generation import nebula_contents
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics import stellar_evolution as evolution
from planetgen.util import draw
from planetgen.web.maps.starmap import uncharted_scene_data
from tests.test_galactic_motion import EDGE_PC, SCATTER_RADIUS_PC, _seed_scatter_grid
from tests.test_uncharted_phenomenon import db_app  # noqa: F401
from planetgen.galaxy.spherical import sector_address_at


def _params(k, age_gy=0.002):
    mass = float(np.linspace(8.0, 40.0, 400)[k])
    return evolution.star_params(mass, age_gy, evolution.evolve_star(mass, age_gy))


def _row(k, params=None):
    params = params or _params(k)
    return {"star_type": params["type"], "yerkes_class": params["yerkes_class"], "mass_kg": params["mass_kg"],
            "radius_km": params["radius_km"], "temperature_k": params["temperature_k"],
            "luminosity_w": params["luminosity_w"], "age_gy": params["age_gy"], "lifespan_gy": params["lifespan_gy"],
            "initial_mass_sol": params["initial_mass_sol"], "phase_end_age_gy": params["phase_end_age_gy"]}


def _built_star(row, k):
    from planetgen.generation import bright_stars as brightStarsGen
    config = SystemConfig()
    config.POPULATION = "young"
    with draw.bound(1000 + k):
        return StarSystem(system_config=config, primary_star_params=brightStarsGen.star_params(row)).stars[0]


def test_a_pre_placed_star_rolls_the_nebula_its_sector_builds():
    """Every class a hot star can light comes out of the row the same as out of the built star."""
    classes = set()
    for k in range(400):
        params = _params(k)
        if evolution.main_sequence_lifetime_gy(params["initial_mass_sol"]) <= params["age_gy"]:
            continue
        row = _row(k, params)
        built = nebula_contents.hosted_nebula_of(_built_star(row, k))
        found = nebula_contents.row_nebula(row)
        if built is None:
            assert found is None
            continue
        assert found is not None
        assert (found.nebula_class, found.hosts, found.radius_ly, found.age_myr, str(found.name)) == (
            built.nebula_class, built.hosts, built.radius_ly, built.age_myr, str(built.name))
        classes.add(found.nebula_class)
    assert classes, "no star of 400 lit a nebula"


def test_the_chance_factor_of_the_cell_scales_the_roll_as_the_sector_does():
    lit = [k for k in range(400) if nebula_contents.row_nebula(_row(k)) is not None]
    assert lit
    assert all(nebula_contents.row_nebula(_row(k), 0.0) is None for k in lit)
    assert all(nebula_contents.row_nebula(_row(k), 1.0) is not None for k in lit)


def _store_hot_rows(mysql_config, count=40, built=()):
    """`count` hot-star rows in the cell of the grid's origin shell; returns the cell's address and the row ids."""
    _seed_scatter_grid(mysql_config)
    point = (SCATTER_RADIUS_PC, 0.0, 0.0)
    address = sector_address_at(point, EDGE_PC)
    conn = store.get_connection(mysql_config)
    ids = []
    try:
        with conn:
            for k in range(count):
                row = _row(k * 5)
                ids.append(conn.execute(
                    "INSERT INTO bright_stars (shell_index, band_index, wedge_index, position_x_mpc, position_y_mpc,"
                    " position_z_mpc, population, star_type, yerkes_class, mass_kg, radius_km, temperature_k,"
                    " luminosity_w, age_gy, initial_mass_sol, seed) VALUES (?, ?, ?, ?, ?, ?, 'young', ?, ?, ?, ?, ?, ?,"
                    " ?, ?, ?)",
                    (*address, point[0] * 1000.0 + k, point[1] * 1000.0, point[2] * 1000.0, row["star_type"],
                     row["yerkes_class"], row["mass_kg"], row["radius_km"], row["temperature_k"], row["luminosity_w"],
                     row["age_gy"], row["initial_mass_sol"], k + 1)).lastrowid)
    finally:
        conn.close()
    return address, ids


def _read(config, function, *args):
    conn = store.get_connection(config)
    try:
        return function(conn, *args)
    finally:
        conn.close()


@pytest.fixture
def gas_rich(monkeypatch):
    """The cell's `"nebula"` factor is 1 (the toy grid's cell lies where the model has next to no gas)."""
    monkeypatch.setattr(hostednebula, "sector_profile", lambda address, skeleton: types.SimpleNamespace(
        kind_factors={"nebula": 1.0}))


def _expected(config, address, row_ids, factor=1.0):
    conn = store.get_connection(config)
    try:
        rows = conn.execute("SELECT * FROM bright_stars WHERE id IN (%s) ORDER BY id" % ",".join("?" * len(row_ids)),
                            row_ids).fetchall()
        return {row["id"]: nebula_contents.row_nebula(row, factor) for row in rows}
    finally:
        conn.close()


def test_a_star_page_detail_carries_the_nebula_of_its_row_and_none_once_built(gas_rich, mysql_config):
    address, ids = _store_hot_rows(mysql_config)
    expected = _expected(mysql_config, address, ids)
    lit = [i for i in ids if expected[i] is not None]
    assert lit and len(lit) < len(ids)
    for row_id in ids:
        detail = _read(mysql_config, query.uncharted_system_detail, row_id)
        if expected[row_id] is None:
            assert detail["hosted_nebula"] is None
        else:
            assert detail["hosted_nebula"] == hostednebula.summary(expected[row_id])
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            conn.execute("SET FOREIGN_KEY_CHECKS = 0")
            conn.execute("UPDATE bright_stars SET star_system_id = 1 WHERE id = ?", (lit[0],))
    finally:
        conn.close()
    assert _read(mysql_config, query.uncharted_system_detail, lit[0])["hosted_nebula"] is None


def test_an_uncharted_cell_lists_the_nebulae_of_its_waiting_stars_and_the_scene_draws_them(gas_rich, mysql_config):
    address, ids = _store_hot_rows(mysql_config)
    expected = _expected(mysql_config, address, ids)
    contents = _read(mysql_config, query.uncharted_sector_contents, *address)
    hosted = contents["hosted_nebulae"]
    assert [entry["star_id"] for entry in hosted] == [i for i in ids if expected[i] is not None]
    for entry in hosted:
        assert entry["name"] == str(expected[entry["star_id"]].name)
        assert entry["radius_ly"] == expected[entry["star_id"]].radius_ly
    center = tuple(sum(entry[axis] for entry in hosted) / len(hosted) for axis in "xyz")
    scene = uncharted_scene_data(address, "T", center, EDGE_PC, contents,
                                 link_url=lambda endpoint, **kw: f"/{endpoint}/{kw.get('bright_star_id') or kw.get('scatter_id')}")
    clouds = [c for c in scene["clouds"] if c["key"].startswith("nebula:h")]
    assert len(clouds) == len(hosted)
    assert all(c["uncharted"] and "nebulaId" not in c and c["href"].startswith("/uncharted_star/") for c in clouds)
    assert all("(uncharted)" in c["name"] for c in clouds)


def test_the_star_page_names_the_nebula_it_lights(gas_rich, mysql_config, db_app):
    address, ids = _store_hot_rows(mysql_config)
    expected = _expected(mysql_config, address, ids)
    lit = next(i for i in ids if expected[i] is not None)
    dark = next(i for i in ids if expected[i] is None)
    page = db_app.test_client().get(f"/systems/uncharted/{lit}").get_data(as_text=True)
    assert "Lights a nebula" in page and str(expected[lit].name) in page
    assert "Lights a nebula" not in db_app.test_client().get(f"/systems/uncharted/{dark}").get_data(as_text=True)


def test_a_cell_with_no_hot_star_lists_no_nebula(mysql_config):
    _seed_scatter_grid(mysql_config)
    assert _read(mysql_config, query.uncharted_sector_contents, 3, 1, 1)["hosted_nebulae"] == []


def test_a_cell_with_next_to_no_gas_lights_nothing(mysql_config):
    """The real profile of the toy grid's cell: the sector's own `"nebula"` factor decides, as it does when it builds."""
    address, ids = _store_hot_rows(mysql_config)
    assert _expected(mysql_config, address, ids, factor=1.0) and any(
        nebula is not None for nebula in _expected(mysql_config, address, ids).values())
    assert _read(mysql_config, query.uncharted_sector_contents, *address)["hosted_nebulae"] == []
    assert all(_read(mysql_config, query.uncharted_system_detail, i)["hosted_nebula"] is None for i in ids)
