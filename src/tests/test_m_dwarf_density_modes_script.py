# tests/test_m_dwarf_density_modes_script.py

"""
GEN.1839: `scripts/m_dwarf_density_modes.py` puts the generator's M-dwarf planets on the density scale of [R2774] (the Earth-like planet of the same mass),
splits them at 0.35 and 0.65 and compares the groups with the published modes. The scale, the body filter, the split and the gap chance are held on small
fixed inputs; the draw itself is `class_census.draw_one`, held by its own test.
"""

import importlib.util
import pathlib

import pytest

from planetgen.physics import bulk_structure

SCRIPT = pathlib.Path(__file__).resolve().parents[2] / "scripts" / "m_dwarf_density_modes.py"


@pytest.fixture(scope="module")
def modes():
    spec = importlib.util.spec_from_file_location("m_dwarf_density_modes", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _body(mass_earth, radius_earth, moon=False, kind="t", host=0):
    return {"is_moon": moon, "body_type": kind, "host": host, "mass": mass_earth * 5.972e24, "radius": radius_earth * 6371.0, "period": 0.05,
            "planet_type": "terrestrial", "water_mass_fraction": 0.01, "envelope_mass_fraction": 0.0}


def _record(*bodies):
    return {"seed": 1, "stars": [{"mass": 0.4 * 1.989e30, "radius": 0.4 * 6.957e5, "temperature": 3400.0}], "bodies": list(bodies)}


@pytest.mark.parametrize("mass", [0.5, 1.0, 4.0, 9.0])
def test_an_earth_like_planet_is_one_at_any_mass(modes, mass):
    radius = bulk_structure.radius_fit_earth(mass, modes.EARTH_LIKE_IRON, 0.0, 1.0)
    assert modes.scaled_density(mass, radius) == pytest.approx(1.0)


def test_a_larger_radius_is_less_dense_by_the_cube(modes):
    radius = bulk_structure.radius_fit_earth(2.0, modes.EARTH_LIKE_IRON, 0.0, 1.0)
    assert modes.scaled_density(2.0, 1.26 * radius) == pytest.approx(0.5, rel=0.01)


def test_only_small_planets_of_a_star_are_kept(modes):
    rows = modes.planets([_record(_body(1.0, 1.0), _body(1.0, 1.0, moon=True), _body(1.0, 1.0, kind="a"), _body(1.0, 4.5), _body(0.1, 0.6),
                                  _body(1.0, 1.0, host=-1), _body(3.0, 1.6))])
    assert [round(row["mass_earth"], 6) for row in rows] == [1.0, 3.0]
    assert rows[0]["host_mass_solar"] == pytest.approx(0.4)
    assert rows[0]["period_d"] == pytest.approx(0.05 * 365.25)


def test_a_failed_system_is_skipped(modes):
    assert modes.planets([{"seed": 2, "error": "boom"}]) == []


def test_the_groups_split_at_the_published_gap(modes):
    rows = [{"density": value, "radius_earth": 1.0, "mass_earth": 1.0} for value in (0.1, 0.34, 0.35, 0.64, 0.65, 1.2)]
    found = modes.groups(rows)
    assert [len(found[name]) for name in ("puffy", "water", "rocky")] == [2, 2, 2]
    assert modes.summary([]) is None


def test_the_chance_of_an_empty_gap_falls_with_the_sample(modes):
    rows = [{"density": 0.65 + 0.2 * (i - 50) / 50.0} for i in range(101)]
    assert 0.0 < modes.gap_chance(rows, 34) < modes.gap_chance(rows, 10) < 1.0


def test_the_published_modes_are_the_three_of_the_paper(modes):
    assert {name: value[0] for name, value in modes.PUBLISHED.items()} == {"rocky": 0.94, "water": 0.47, "puffy": 0.24}
