"""GEN.1221: the class-by-domain matrix and the quantities a body gives it."""

import types

import pytest

from planetgen.physics import class_domains as domains, class_fit, planet_classes


def _values(**fields):
    base = {"surface_temperature_k": 288.0, "atmospheric_pressure_pa": 101325.0, "hydrosphere": "surface ocean", "relative_humidity": 0.6, "groundwater": "none"}
    for number, mean in enumerate((250.0, 270.0, 285.0, 295.0, 285.0, 270.0, 250.0), start=1):
        base[f"band{number}_mean_k"], base[f"band{number}_swing_k"] = mean, 10.0
    base.update(fields)
    return base


def test_regional_day_and_night_temperatures_come_from_the_seven_bands():
    found = domains.derived(_values())
    assert found["hottest_region_k"] == 295.0 and found["coldest_region_k"] == 250.0
    assert found["day_temperature_k"] == 300.0 and found["night_temperature_k"] == 245.0


def test_a_body_without_bands_has_no_regional_temperature():
    values = _values()
    values["band3_mean_k"] = None
    found = domains.derived(values)
    assert found["hottest_region_k"] is None and found["wet_bulb_c"] is None


def test_the_wet_bulb_is_below_the_hottest_day_and_follows_the_humidity():
    dry = domains.derived(_values(relative_humidity=0.1))["wet_bulb_c"]
    wet = domains.derived(_values(relative_humidity=0.95))["wet_bulb_c"]
    assert dry < wet < 300.0 - 273.15 + 0.5


def test_groundwater_and_the_phase_gates_are_one_or_zero():
    found = domains.derived(_values(groundwater="reachable"))
    assert found["groundwater_present"] == 1.0 and domains.derived(_values())["groundwater_present"] == 0.0
    assert found["phase_water"] == 1.0 and found["phase_co2"] == 0.0 and found["phase_methane_ethane"] == 0.0


def test_every_quantity_of_the_fit_table_has_a_domain_and_every_class_a_cell_in_each():
    for cls in class_fit.CLASSES:
        cells = domains.matrix(cls)
        assert set(cells) == set(domains.ALL_DOMAINS)
        assert sum(len(conditions) for conditions in cells.values()) == len(cls.when)


def test_the_oxygen_water_family_has_the_ranges_of_the_design():
    by_letter = {cls.letter: cls for cls in class_fit.CLASSES if cls.kind == "t"}
    for letter in "MHOPRZUY":
        assert letter in by_letter, letter
    thermal = {letter: [c for c in domains.matrix(by_letter[letter])["thermal"]] for letter in "MHOPR"}
    assert {"quantity": "temperature_k", "min": 275.0} in thermal["M"] and {"quantity": "temperature_k", "lt": 305.0} in thermal["M"]
    assert {"quantity": "temperature_k", "min": 312.0} in thermal["R"]
    assert {"quantity": "wet_bulb_c", "min": 31.0} in domains.matrix(by_letter["R"])["wet_bulb"]


def test_an_earth_like_body_scores_one_in_every_domain_of_m():
    values = _values()
    values.update(body_type="t", mass_kg=5.97e24, gravity_g=1.0, ocean_fraction=0.7, water_mass_fraction=2e-4)
    values["gas_partials_pa"] = {"N2": 78000.0, "O2": 21000.0, "Ar": 900.0, "CO2": 40.0}
    found = class_fit._found(values)
    scores = domains.domain_scores(next(c for c in class_fit.CLASSES if c.letter == "M" and c.kind == "t"), found, class_fit.membership)
    assert all(score == pytest.approx(1.0, abs=0.01) for score in scores.values()), scores
    assert class_fit.best_class(values)[0] == "M"


def test_overall_is_the_weighted_geometric_mean_and_zero_when_a_domain_is_zero():
    weights = {domain: 1.0 for domain in domains.DOMAINS}
    scores = {domain: 0.5 for domain in domains.ALL_DOMAINS}
    assert domains.overall(scores, weights) == pytest.approx(0.5)
    scores["water"] = 0.0
    assert domains.overall(scores, weights) == 0.0
