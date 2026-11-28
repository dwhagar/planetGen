"""scripts/gen1554_zone_edge_sweep.py (GEN.1554): the sweep lists the primary star's planets through its host, skips belts, and its report reads what it kept."""

import importlib.util
import pathlib
from types import SimpleNamespace

_PATH = pathlib.Path(__file__).resolve().parents[2] / "scripts" / "gen1554_zone_edge_sweep.py"
_SPEC = importlib.util.spec_from_file_location("gen1554_zone_edge_sweep", _PATH)
sweep = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(sweep)

_EARTH = sweep.EARTH_MASS_KG


def _planet(mass_earth, distance, **extra):
    fields = dict(body_type="t", mass=mass_earth * _EARTH, distance=distance, moons=[], gas_partials_pa={"CO2": 4.0e4, "N2": 8.0e4},
                  surface_temperature=288.0, atmospheric_pressure=1.2e5, bond_albedo=0.3, greenhouse_factor=33.0,
                  ocean_fraction=0.7, planet_class="M", eccentricity=0.02)
    fields.update(extra)
    return SimpleNamespace(**fields)


def _system(orbiters):
    star = SimpleNamespace(habitable_zone=(0.95, 1.68), luminosity=sweep.SOLAR_LUMINOSITY_W, temperature=5772, age=4.6)
    return SimpleNamespace(primary_star=star, hosts=[SimpleNamespace(star=star, role="primary", orbiters=orbiters)])


def test_the_planets_come_from_the_primary_host_and_belts_are_left_out():
    belt = SimpleNamespace(body_type="a", mass=None, distance=1.0, moons=[])
    rows = sweep._bodies_of(_system([belt, _planet(1.0, 1.0)]), "G2V", 7)
    assert len(rows) == 1
    row = rows[0]
    assert (row["star"], row["index"], row["moon"], row["cls"]) == ("G2V", 7, False, "M")
    assert abs(row["flux"] - 1.0) < 1e-9
    assert 0.0 < row["position"] < 1.0


def test_a_system_with_no_planets_gives_no_rows_and_no_error():
    assert sweep._bodies_of(_system([]), "M5V", 0) == []
    assert sweep._bodies_of(SimpleNamespace(primary_star=_system([]).primary_star, hosts=[]), "M5V", 1) == []


def test_a_body_outside_the_mass_or_flux_window_is_not_kept():
    light = _planet(0.1, 1.0)
    heavy = _planet(9.0, 1.0)
    far = _planet(1.0, 40.0)
    assert sweep._bodies_of(_system([light, heavy, far]), "G2V", 0) == []


def test_a_moon_of_the_right_size_is_kept_beside_its_planet():
    moon = _planet(0.5, 1.0)
    giant = _planet(300.0, 1.0, moons=[moon])
    rows = sweep._bodies_of(_system([giant]), "G2V", 0)
    assert [r["moon"] for r in rows] == [True]


def test_the_report_has_a_flux_table_and_a_position_table_for_each_star():
    rows = sweep._bodies_of(_system([_planet(1.0, 1.0), _planet(2.0, 1.4)]), "G2V", 0)
    text = sweep.report_text(rows)
    assert "## G2V (2 bodies)" in text
    assert "By incident flux" in text and "By position in the zone" in text
    assert "| 1 | 1 | 288 | 1 of 1 |" in text
