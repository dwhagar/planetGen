# tests/test_sweep_atmosphere.py

"""
TEST.209: one consistency sweep over every planet and moon class. Planets of every class and zone round three host stars, and the
planets and moons of generated systems, are checked by the rules of `tests/audit/sweep_atmosphere.py`: the partial pressures add up to the
pressure; the molar mass, density and scale height follow from them; the temperature is the formula on the stored albedo and greenhouse
quantity; the ozone follows the oxygen and the star; the aerosol depths follow their rules; every gas in the text is in the table and the
mixture. The last test breaks each rule in turn and checks that the sweep reports it.

The first half tests the sweep's own frame on synthetic bodies.
"""

import types

import pytest

from tests.audit import sweep
from tests.audit import sweep_atmosphere as rules


def _entry(label="x", **attributes):
    return sweep.Entry(types.SimpleNamespace(**attributes), None, None, label)


# --- The frame ---------------------------------------------------------------------------------------------------------------


def test_a_rule_that_holds_reports_nothing_and_one_that_fails_names_the_entry():
    rule = sweep.Rule("positive", lambda entry: [] if entry.body.value > 0 else [f"value {entry.body.value}"],
                      lambda entry: sweep.broken_copy(entry.body, value=-1.0))
    found = sweep.run_rules((rule,), (_entry("good", value=1.0), _entry("bad", value=-2.0)))
    assert found == {"positive": [("bad", "value -2.0")]}


def test_a_rule_whose_breaker_is_caught_can_fail_and_one_that_is_not_is_reported_false():
    catches = sweep.Rule("catches", lambda entry: [] if entry.body.value > 0 else ["bad"],
                         lambda entry: sweep.broken_copy(entry.body, value=-1.0))
    misses = sweep.Rule("misses", lambda entry: [], lambda entry: sweep.broken_copy(entry.body, value=-1.0))
    skips = sweep.Rule("skips", lambda entry: [], lambda entry: None)
    result = sweep.prove_rules_can_fail((catches, misses, skips), (_entry(value=1.0),))
    assert result == {"catches": True, "misses": False, "skips": None}


def test_a_broken_copy_never_touches_the_original():
    body = types.SimpleNamespace(value=1.0, gases={"N2": 1.0})
    clone = sweep.broken_copy(body, value=2.0)
    assert body.value == 1.0 and clone.value == 2.0 and clone.gases is body.gases


def test_every_rule_has_a_name_a_check_and_a_breaker_and_the_names_are_unique():
    assert len({rule.name for rule in rules.RULES}) == len(rules.RULES) >= 11
    assert all(callable(rule.check) and callable(rule.breaker) for rule in rules.RULES)


# --- The sweep -------------------------------------------------------------------------------------------------------------------


@pytest.fixture(scope="module")
def everything():
    return sweep.everything()


@pytest.mark.slow
def test_the_sweep_holds_planets_of_many_classes_in_every_zone_and_some_moons(everything):
    built = {(entry.body.planet_class, entry.body.zone) for entry in everything if not entry.is_moon}
    assert {zone for _letter, zone in built} >= set(sweep.ZONES)
    assert len({letter for letter, _zone in built}) >= 6      # GEN.646: classes fall where the conditions put them
    assert sum(1 for entry in everything if entry.is_moon) >= 20
    assert {entry.star.type[:3] for entry in everything} >= {"M2V", "G2V", "F5V"}


@pytest.mark.slow
def test_the_sweep_repeats(everything):
    again = sweep.entries.__wrapped__()
    first = [(entry.label, entry.body.surface_temperature, entry.body.atmospheric_pressure) for entry in everything[:len(again)]]
    assert first == [(entry.label, entry.body.surface_temperature, entry.body.atmospheric_pressure) for entry in again]


@pytest.mark.slow
def test_every_body_of_the_sweep_obeys_every_rule(everything):
    found = sweep.run_rules(rules.RULES, everything)
    failures = {name: violations[:5] for name, violations in found.items() if violations}
    assert not failures, "\n".join(f"{name}: {len(found[name])} bodies, first {violations}" for name, violations in failures.items())


@pytest.mark.slow
def test_every_rule_is_exercised_by_air_that_the_sweep_holds(everything):
    airy = [entry for entry in everything if rules._air(entry.body)]
    assert len(airy) >= 80      # enough air-holding bodies to exercise every rule (97 of the sweep on 11-01; the draws move it a little)
    assert any(entry.body.tau_cloud for entry in airy if entry.body.tau_cloud)
    assert any(not rules._air(entry.body) for entry in everything)


def _oxygen_entry(column_du, p_o2_pa=21200.0, teff=5800.0):
    star = types.SimpleNamespace(temperature=teff)
    body = types.SimpleNamespace(gas_partials_pa={"N2": 79000.0, "O2": p_o2_pa}, gravity=1.0, ozone_column_du=column_du)
    return sweep.Entry(body, None, star, "oxygen")


def test_the_temperature_rule_stops_a_runaway_greenhouse_at_the_molten_surface_limit(monkeypatch):
    """Stand-in CC7: the air solve holds the surface at `MAX_SURFACE_TEMPERATURE_K`, so a hotter greenhouse still gives that, not the Stefan-Boltzmann value."""
    from planetgen.physics import surface_air

    monkeypatch.setattr(rules, "_starlight_w_m2", lambda entry: 1.0e5)
    entry = _entry("hot", bond_albedo=0.1, greenhouse_factor=100.0, surface_temperature=surface_air.MAX_SURFACE_TEMPERATURE_K, gas_partials_pa={"N2": 1.0e5})
    assert rules._check_temperature(entry) == []
    entry.body.surface_temperature = surface_air.MAX_SURFACE_TEMPERATURE_K * 1.1
    assert rules._check_temperature(entry)


def test_the_ozone_rule_lets_a_body_that_flares_hold_the_flag_keep_a_tenth_of_the_least_column():
    """GEN.808: a world whose star's flares hold the ozone flag keeps `OZONE_LOSS_REMAINING` of its column, so the least the rule allows falls by that."""
    from planetgen.physics import ozone, radiation

    median = ozone.median_column_du(5800.0, 1.0)
    low, _high = ozone.structure_range(ozone.o2_level_pal(21200.0, 1.0))
    column = median * low * radiation.OZONE_LOSS_REMAINING * 1.05
    plain = _oxygen_entry(column)
    flagged = _oxygen_entry(column)
    flagged.body.ozone_loss_flag = True
    assert rules._check_ozone(plain) and rules._check_ozone(flagged) == []
    flagged.body.ozone_column_du = median * low * radiation.OZONE_LOSS_REMAINING * 0.9
    assert rules._check_ozone(flagged)


def test_the_ozone_rule_holds_for_earth_air_and_fails_outside_the_fits_spread():
    """The sweep's bodies hold no free oxygen since the air solve of 4c (GEN.855), so the rule's oxygen branch is tried on a body
    made for it: Earth's oxygen under the Sun gives a column in the fit's spread, 100 times that and none at all do not."""
    from planetgen.physics import ozone

    median = ozone.median_column_du(5800.0, 1.0)
    low, high = ozone.structure_range(ozone.o2_level_pal(21200.0, 1.0))      # the spread of the models about the fit at this oxygen level
    assert rules._check_ozone(_oxygen_entry(median * 0.5 * (low + high))) == []
    assert rules._check_ozone(_oxygen_entry(median * high * 1.001))
    assert rules._check_ozone(_oxygen_entry(median * low * 0.999))
    assert rules._check_ozone(_oxygen_entry(1.0, p_o2_pa=0.0))      # ozone with no oxygen
    assert rules._check_ozone(_oxygen_entry(0.0)) != []              # no ozone with Earth's oxygen
    low, high = ozone.abiotic_column_range_du(5800.0)
    assert rules._check_ozone(_oxygen_entry(high, p_o2_pa=0.002)) == []      # abiotic oxygen, below 0.01 percent of the present level


@pytest.mark.slow
def test_each_rule_reports_a_body_with_exactly_that_rule_broken(everything):
    result = sweep.prove_rules_can_fail(rules.RULES, everything)
    assert result == {rule.name: True for rule in rules.RULES}, result
