"""GEN.442: the mixing rule and the tapered brittle compaction of small bodies (`planetgen/physics/compaction.py`, its table
`physics/data/compaction_table.json`, built by `scripts/build_compaction_table.py`; docs/design/small-bodies-composition.md
section 4). The Solar System bodies below are test targets only."""

import importlib.util
import json
import math
import os

import pytest

from planetgen.physics import compaction as cp
from planetgen.physics.elements import DataError

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..")


def _builder():
    spec = importlib.util.spec_from_file_location("build_compaction_table", os.path.join(ROOT, "scripts", "build_compaction_table.py"))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _density(mass_kg, rock):
    return mass_kg / (4.0 / 3.0 * math.pi * (cp.compacted_radius_km(mass_kg, rock) * 1e3) ** 3)


# --- the laws ----------------------------------------------------------------------------------------------------

def test_the_brittle_laws_are_bierson_and_nimmo_table_1():
    # the laws are in the table builder (the table cannot be read before it is built); the module keeps the grain densities
    build = _builder()
    assert build.ice_porosity(0.0) == pytest.approx(0.60)
    assert build.ice_porosity(5e6) == pytest.approx(0.60 * math.exp(-0.5))
    assert build.ice_porosity(1e9) == pytest.approx(0.20)
    assert build.silicate_porosity(1e5) == pytest.approx(0.6)
    assert build.silicate_porosity(1e8) == pytest.approx(0.6 * 100.0 ** -0.11)
    assert cp.ICE_GRAIN_KG_M3 == 920.0 and cp.SILICATE_GRAIN_KG_M3 == 3500.0


def test_the_taper_is_one_below_1e19_zero_above_1e21_and_linear_in_log_mass():
    assert cp.taper_weight(1e15) == 1.0 and cp.taper_weight(1e19) == 1.0
    assert cp.taper_weight(1e21) == 0.0 and cp.taper_weight(1e24) == 0.0
    assert cp.taper_weight(1e20) == pytest.approx(0.5)
    assert cp.taper_weight(6.3e20) == pytest.approx((21 - math.log10(6.3e20)) / 2.0)
    with pytest.raises(ValueError):
        cp.taper_weight(0.0)


def test_the_mixing_rule_is_one_over_the_sum_of_fractions_over_densities():
    rho = cp.grain_density_kg_m3({"water ice": 0.3, "forsterite": 0.7})
    assert rho == pytest.approx(1.0 / (0.3 / 917.0 + 0.7 / 3271.0))
    # fractions are scaled, a stated density overrides the table, and the result lies between the grains
    assert cp.grain_density_kg_m3({"water ice": 3.0, "forsterite": 7.0}) == pytest.approx(rho)
    assert cp.grain_density_kg_m3({"complex organics": 1.0}, densities={"complex organics": 2000.0}) == pytest.approx(2000.0)
    assert 917.0 < rho < 3271.0


def test_a_bad_mixture_is_refused():
    for call in (lambda: cp.grain_density_kg_m3({}), lambda: cp.grain_density_kg_m3({"forsterite": -1.0, "iron": 2.0}),
                 lambda: cp.grain_density_kg_m3({"olivine": 1.0})):
        with pytest.raises(ValueError):
            call()
    with pytest.raises(KeyError):
        cp.grain_density_kg_m3({"unobtainium": 1.0})


# --- the table against the reference integration ------------------------------------------------------------------

def test_the_design_notes_bodies_come_out_of_the_table():
    # rock fraction 0.70 (taper.py of the design note): Orcus-sized 6.3e20 kg has a radius near 433 km and 1.85 g/cm3
    assert cp.compacted_radius_km(6.3e20, 0.70) == pytest.approx(433.0, rel=0.02)
    assert _density(6.3e20, 0.70) == pytest.approx(1847.0, rel=0.02)
    assert _density(1e21, 0.70) == pytest.approx(1907.0, rel=0.02)
    assert _density(3e20, 0.70) == pytest.approx(1739.0, rel=0.02)
    assert _density(1e19, 0.70) == pytest.approx(815.0, rel=0.04)
    assert _density(1e20, 0.70) == pytest.approx(1526.0, rel=0.02)


def test_the_table_agrees_with_the_integration_at_a_few_masses():
    build = _builder()
    for mass, rock in ((3e16, 0.55), (4e19, 0.80), (5e20, 0.45)):
        weight = build.taper(mass)
        curve = build.density_curve(rock, weight)
        import numpy as np
        log_p = np.log10(build.PRESSURE_GRID_PA)

        def density_of_pressure(p, curve=curve, log_p=log_p):
            return float(np.interp(math.log10(max(p, 1e-3)), log_p, curve))

        from scipy.optimize import brentq
        log_pc = brentq(lambda x: math.log(build.sphere(10.0 ** x, density_of_pressure)[0]) - math.log(mass), -2.5, 9.6, xtol=1e-6)
        got_mass, radius_m = build.sphere(10.0 ** log_pc, density_of_pressure)
        assert cp.compacted_radius_km(mass, rock) == pytest.approx(radius_m / 1e3, rel=0.02), (mass, rock)
        assert cp.compacted_density_kg_m3(radius_m / 1e3, rock) == pytest.approx(got_mass / (4.0 / 3.0 * math.pi * radius_m ** 3),
                                                                                  rel=0.02), (mass, rock)


def test_the_density_is_continuous_in_mass_across_the_taper_and_never_falls_with_mass():
    for rock in (0.35, 0.55, 0.70, 0.85):
        previous = 0.0
        masses = [10.0 ** (12.0 + 0.05 * i) for i in range(181)]
        for mass in masses:
            value = _density(min(mass, 1e21), rock)
            assert value >= previous * 0.995, (rock, mass)
            previous = value
        for edge in (1e19, 1e21):
            step = abs(_density(edge * 1.1 if edge < 1e21 else edge, rock) - _density(edge * 0.9, rock)) / _density(edge * 0.9, rock)
            assert step < 0.2, (rock, edge, step)


def test_rockier_bodies_are_denser_and_icier_ones_larger():
    for mass in (1e15, 1e19, 5e20):
        densities = [_density(mass, rock) for rock in (0.3, 0.5, 0.7, 0.9)]
        assert densities == sorted(densities)
        radii = [cp.compacted_radius_km(mass, rock) for rock in (0.3, 0.5, 0.7, 0.9)]
        assert radii == sorted(radii, reverse=True)


def test_radius_and_mass_round_trip():
    for mass, rock in ((1e14, 0.6), (2e18, 0.4), (7e20, 0.75)):
        radius = cp.compacted_radius_km(mass, rock)
        assert cp.compacted_density_kg_m3(radius, rock) * 4.0 / 3.0 * math.pi * (radius * 1e3) ** 3 == pytest.approx(mass, rel=0.01)


def test_outside_the_table_is_an_error_and_never_a_guess():
    with pytest.raises(cp.OutOfRange):
        cp.compacted_radius_km(1e22, 0.7)
    with pytest.raises(cp.OutOfRange):
        cp.compacted_density_kg_m3(cp.max_radius_km(0.7) * 1.05, 0.7)
    with pytest.raises(ValueError):
        cp.compacted_density_kg_m3(100.0, 0.1)
    with pytest.raises(ValueError):
        cp.compacted_density_kg_m3(-1.0, 0.5)
    assert cp.compacted_density_kg_m3(0.01, 0.7) == pytest.approx(cp.compacted_density_kg_m3(0.7, 0.7), rel=0.02)


# --- the data file -------------------------------------------------------------------------------------------------

def test_the_data_file_states_its_sources_and_is_refused_without_them(monkeypatch, tmp_path):
    with open(cp.DATA_FILE, encoding="utf-8") as handle:
        data = json.load(handle)
    assert set(data["sources"]) == {"bierson-nimmo-2019", "stand-in-eos", "taper"}
    broken = dict(data, sources={})
    path = tmp_path / "compaction_table.json"
    path.write_text(json.dumps(broken), encoding="utf-8")
    monkeypatch.setattr(cp, "DATA_FILE", str(path))
    with pytest.raises(DataError):
        cp._load()
    broken = dict(data, taper_end_kg=1e22)
    path.write_text(json.dumps(broken), encoding="utf-8")
    with pytest.raises(DataError):
        cp._load()
