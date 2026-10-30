# tests/test_shoreline_caps.py

"""GEN.483: the pressure caps of row d6 are the envelope of the Solar System pressures by depth below the Berta-Thompson shoreline (`physics/shoreline_caps.py`)."""

import pytest

from planetgen.physics import shoreline_caps as sc
from planetgen.generation.possibility_table import Subject
from planetgen.physics import atmosphere_recipes as ar
from planetgen.physics.atmosphere_recipes import RECIPES


def body(name):
    return next(b for b in sc.BODIES if b.name == name)


def test_the_shoreline_parameters_are_the_published_ones():
    assert (sc.LOG_F0, sc.P_V, sc.Q_L) == (2.62, 5.9, 1.17)
    # the recipe's own expression of the same line
    with open(ar.DATA_FILE, encoding="utf-8") as handle:
        assert "2.62 + 5.9*log10(v_esc_kms/11.186) + 1.17*log10(luminosity_solar)" in handle.read()


def test_the_earth_sits_2_6_dex_below_the_line_and_a_body_on_it_is_at_zero():
    assert sc.depth(body("Earth")) == pytest.approx(2.62)
    assert sc.delta(1.0, 11.186) == pytest.approx(-2.62)
    # a flux of 10^2.62 is at 1 / sqrt(10^2.62) AU
    assert abs(sc.delta(10 ** (-2.62 / 2), 11.186)) < 0.01


def test_the_solar_system_bodies_fall_where_the_paper_places_them():
    depths = {b.name: sc.depth(b) for b in sc.BODIES}
    assert depths["Pluto"] == pytest.approx(0.11, abs=0.02)
    assert depths["Triton"] == pytest.approx(0.35, abs=0.02)
    assert depths["Titan"] == pytest.approx(0.88, abs=0.02)
    assert depths["Mars"] == pytest.approx(0.94, abs=0.02)
    assert depths["Venus"] == pytest.approx(2.14, abs=0.02)
    assert all(d > 0 for d in depths.values()), "every one is on the atmosphere side"


def test_the_envelope_keeps_only_record_pressures_and_leaves_mars_and_the_earth_under_it():
    names = []
    for depth, level in sc.envelope():
        names.append(next(b.name for b in sc.BODIES if sc.depth(b) == depth))
    assert names == ["Pluto", "Triton", "Titan", "Venus"]
    levels = [level for _, level in sc.envelope()]
    assert levels == sorted(levels)


def test_the_cap_rises_with_depth_and_is_absent_below_the_deepest_record():
    caps = [sc.pressure_cap_bar(d) for d in (0.0, 0.2, 0.5, 1.0, 1.5, 2.0)]
    assert caps == sorted(caps)
    assert sc.pressure_cap_bar(0.0) == pytest.approx(1.3e-5)
    assert sc.pressure_cap_bar(2.2) is None
    assert sc.pressure_cap_bar(sc.depth(body("Titan"))) == pytest.approx(1.467)


def test_every_measured_body_is_under_the_cap_of_its_own_depth():
    for b in sc.BODIES:
        assert b.pressure_bar <= sc.pressure_cap_bar(sc.depth(b)) * (1 + 1e-9) if sc.pressure_cap_bar(sc.depth(b)) else True, b.name


def test_the_recipe_takes_the_caps_the_window_edges_give():
    """The window edges of rows d2 to d5 are -0.3, -0.8 and -1.4 dex (a design); the d6 table is the envelope read at their depths."""
    caps = sc.window_caps((0.3, 0.8, 1.4))
    assert caps == [1.4e-05, 0.25, 8.0, None]
    for row, cap in zip(("d2", "d3", "d4", "d5"), caps):
        assert RECIPES.derive("d6", Subject("terrestrial", {"drawn.d": row})) == {"p_max_bar": cap}
    assert RECIPES.derive("d6", Subject("terrestrial", {"drawn.d": "d1"})) == {"p_max_bar": 0.0}
