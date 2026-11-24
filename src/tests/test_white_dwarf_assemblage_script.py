# tests/test_white_dwarf_assemblage_script.py

"""
GEN.1841: `scripts/white_dwarf_assemblage.py` puts the mantles of simulated rocky bodies on the mineral norm and counts the exotic ones (quartz or periclase) the way the polluted
white dwarfs are read. The oxide readings, the exotic test, the interval of the published count and the draw (deterministic, a norm that adds to 1) are held on small inputs.
"""

import importlib.util
import pathlib

import pytest

SCRIPT = pathlib.Path(__file__).resolve().parents[2] / "scripts" / "white_dwarf_assemblage.py"


@pytest.fixture(scope="module")
def wd():
    spec = importlib.util.spec_from_file_location("white_dwarf_assemblage", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def drawn(wd):
    return wd.bodies(40, 4242)


def test_exotic_needs_quartz_or_periclase_at_the_threshold(wd):
    assert wd.exotic({"quartz": 0.1, "olivine": 0.9}, 0.10)
    assert wd.exotic({"periclase": 0.25}, 0.20)
    assert not wd.exotic({"quartz": 0.04, "periclase": 0.049}, 0.05)
    assert not wd.exotic({}, 0.05)


def test_the_published_counts_are_the_stated_ones(wd):
    assert wd.PUBLISHED["putirka-xu"] == (12, 23)
    assert wd.PUBLISHED["trierweiler-raw-not-chondritic"] == (31 - 15, 31)
    assert wd.PUBLISHED["trierweiler-steady-state-not-chondritic"] == (31 - 21, 31)


def test_the_wilson_interval_of_twelve_in_twenty_three(wd):
    low, high = wd.wilson(12, 23)
    assert low == pytest.approx(0.33, abs=0.01) and high == pytest.approx(0.71, abs=0.01)
    assert low < 12 / 23 < high


def test_every_reading_gives_a_norm_that_adds_to_one(wd, drawn):
    for body in drawn:
        assert set(body["shares"]) == set(wd.READINGS)
        for shares in body["shares"].values():
            assert sum(shares.values()) == pytest.approx(1.0, abs=1e-9)


def test_the_draw_is_deterministic_and_the_stars_are_stellar(wd, drawn):
    again = wd.bodies(40, 4242)
    assert [b["mg_si"] for b in again] == [b["mg_si"] for b in drawn]
    assert all(0.6 < b["mg_si"] < 2.0 and 0.03 < b["ca_mg"] < 0.12 and 0.4 < b["fe_mg"] < 1.5 for b in drawn)


def test_a_mantle_with_the_iron_in_the_core_is_a_peridotite_and_with_all_the_iron_is_not(wd, drawn):
    built = sum(wd.exotic(b["shares"]["four"], 0.10) for b in drawn)
    all_iron = sum(wd.exotic(b["shares"]["all-iron"], 0.10) for b in drawn)
    assert built <= 1 and all_iron > 5


def test_the_quantiles_are_ordered(wd):
    low, middle, high = wd.quantiles([5, 1, 3, 2, 4, 9, 8, 7, 6, 10])
    assert low < middle < high
