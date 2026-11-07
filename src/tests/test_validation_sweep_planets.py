"""
The statistics of the validation sweep of planets, atmospheres, radiation and heat (GEN.753, `scripts/validation_sweep_planets.py`) on a population small enough to
work out by hand. The sweep itself (hundreds of generated systems) is run by the script, and its numbers are in `reports/gen753-planets-atmospheres-radiation-heat-2026-11-07.md`.
"""

import pathlib
import sys

import pytest

ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))

import validation_sweep_planets as sweep  # noqa: E402

STAR = {"mass": 2.0e30, "radius": 7.0e5, "temperature": 5800.0, "age": 4.5, "luminosity": 3.828e26}


def _body(distance, mass_earth, radius_earth=1.0, **extra):
    row = {name: None for name in sweep.FIELDS}
    row.update(body_type="t", planet_class="E", mass=mass_earth * sweep.EARTH_MASS_KG, radius=radius_earth * sweep.EARTH_RADIUS_KM, distance=distance,
             period=distance ** 1.5, eccentricity=0.02, atmospheric_pressure=0.0, system=0, planet_index=0, planets_in_system=3, is_moon=False, star=STAR)
    row.update(extra)
    return row


def test_the_spacing_is_in_mutual_hill_radii_of_neighbours_whatever_the_order():
    inner, outer = _body(0.1, 1.0), _body(0.15, 1.0)
    hill = ((2.0 * sweep.EARTH_MASS_KG) / (3.0 * STAR["mass"])) ** (1.0 / 3.0) * 0.125
    assert sweep.mutual_hill_spacing([outer, inner]) == pytest.approx([0.05 / hill])
    assert sweep.mutual_hill_spacing([inner]) == []
    assert sweep.mutual_hill_spacing([]) == []


def test_the_summary_counts_the_valley_the_spacing_and_the_dose_of_a_small_population():
    rows = [_body(0.05, 3.0, 1.6), _body(0.08, 8.0, 2.4), _body(0.12, 1.0, 1.0), _body(1.0, 1.0, 1.0, dose_gcr_msv_yr=400.0, planet_index=1)]
    summary = sweep.summarize(rows)
    assert summary["planets"] == 4 and summary["systems"] == 1
    assert summary["radius_histogram"]["n"] == 3                            # the three planets inside 100 days, 1 to 4 Earth radii
    assert summary["radius_histogram"]["share_over_2"] == pytest.approx(1.0 / 3.0)
    assert summary["spacing"]["n"] == 3 and summary["spacing"]["share_over_10"] <= 1.0
    assert summary["airless_dose_at_1_au"]["n"] == 1 and summary["airless_dose_at_1_au"]["median"] == pytest.approx(400.0)
    assert "radius histogram" in sweep.report(summary)


def test_a_population_with_nothing_to_measure_reports_without_failing():
    summary = sweep.summarize([_body(1.0, 1.0, planets_in_system=1)])
    assert summary["spacing"]["n"] == 0 and summary["inner_small_eccentricity"]["n"] == 0
    assert sweep.report(summary)
