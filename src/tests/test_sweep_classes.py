# tests/test_sweep_classes.py

"""
TEST.301: classes from conditions, on the shared body sweep (`tests/audit/sweep.py`). Every planet and moon has a class of the rule file
and of its own kind; `classify` is total and pure on all of them; the stored letter is the one the conditions give (a known anomaly until
the class leaves generation, GEN.646, so the rule reports and does not fail while the tag stands); and a system generated with
`+habitable_world` answers the request (a known anomaly while temperate liquid-water worlds are rare and the request is steered by the letter, GEN.909, GEN.853, GEN.854, GEN.916). Each rule has a broken copy
that the sweep must report.

The first half tests the rules on synthetic bodies.
"""

import json
import pathlib
import types

import pytest

from planetgen.physics import atmosphere, constants, planet_classes as pc
from tests.audit import prevalence as census
from tests.audit import sweep
from tests.audit import sweep_classes as rules

EXPECTATIONS = pathlib.Path(__file__).resolve().parent / "data" / "class_sweep_expectations.json"
REPORT_NAME = "class-sweep-latest.md"


def _body(**changes):
    """Earth's stored conditions as a body object; its class letter is M."""
    fields = dict(
        body_type="t", mass=constants.EARTH_MASS_TO_KG, gravity=1.0, surface_temperature=288.0, atmospheric_pressure=101325.0,
        hydrosphere="surface ocean", ocean_fraction=0.7, water_mass_fraction=2.3e-4, flux_variation=0.0, ocean_ph=None,
        planet_class="M", **{atmosphere.GAS_FIELD: {"O2": 21200.0, "N2": 78800.0, "Ar": 900.0, "CO2": 40.0}})
    fields.update(changes)
    return types.SimpleNamespace(**fields)


def _entry(body=None):
    return sweep.Entry(body or _body(), None, None, "x")


def _rule(name):
    return next(rule for rule in rules.RULES if rule.name == name)


# --- The rules on bodies made for them -------------------------------------------------------------------------------


def test_an_earth_like_body_with_its_conditions_class_obeys_every_rule():
    for rule in rules.RULES:
        assert rule.check(_entry()) == [], rule.name


def test_a_body_without_a_class_or_with_an_unknown_one_or_one_of_the_wrong_kind_is_reported():
    check = _rule("class-assigned").check
    assert check(_entry(_body(planet_class=None))) and check(_entry(_body(planet_class="@")))
    assert check(_entry(_body(planet_class="J")))                 # a giant class on a terrestrial body
    assert check(_entry(_body(body_type="g", planet_class="M")))


def test_a_body_that_cannot_be_classified_is_reported():
    assert _rule("conditions-total").check(_entry(_body(body_type="x")))
    assert _rule("conditions-total").check(_entry(_body(mass=None)))


def test_the_stored_letter_must_be_the_conditions_letter():
    check = _rule("class-from-conditions").check
    assert check(_entry(_body(planet_class="C"))) != []
    dry = _body(atmospheric_pressure=0.0, hydrosphere="dry", ocean_fraction=0.0, planet_class="C", **{atmosphere.GAS_FIELD: {}})
    dry.mass, dry.gravity = 0.15 * constants.EARTH_MASS_TO_KG, 0.3      # Class C from 0.1 Earth masses (GEN.1160)
    assert check(_entry(dry)) == []
    assert check(_entry(_body(planet_class="D"))) and "'D'" in check(_entry(_body(planet_class="D")))[0]


def test_every_rule_has_a_breaker_the_check_reports():
    for rule in rules.RULES:
        broken = rule.breaker(_entry())
        assert broken is not None, rule.name
        assert rule.check(sweep.Entry(broken, None, None, "x")), rule.name


def test_the_conditions_letter_is_the_first_row_of_the_kind_that_holds():
    assert rules._conditions_letter(_body()) == pc.classify(_body()) == "M"
    giant = _body(body_type="g", mass=20.0 * constants.EARTH_MASS_TO_KG)
    assert rules._conditions_letter(giant) == pc.classify(giant) == "I"


def test_the_share_of_met_requests_counts_aim_and_fallback():
    assert rules.request_met_share(["aim", "fallback", "unmet", "unmet"]) == 0.5
    assert rules.request_met_share([]) == 0.0


# --- The sweep ---------------------------------------------------------------------------------------------------------------


@pytest.fixture(scope="module")
def expectations():
    return json.loads(EXPECTATIONS.read_text(encoding="utf-8"))


@pytest.fixture(scope="module")
def everything():
    return sweep.everything()


def test_the_expectations_name_an_item_for_every_anomaly(expectations):
    for name, text in expectations["anomalies"].items():
        assert name in {rule.name for rule in rules.RULES} | {"habitable-request"} and text.split(":")[0].split(" ")[0].startswith("GEN."), name


@pytest.mark.slow
def test_every_body_is_classed_and_classify_is_total_and_pure_on_all_of_them(everything):
    found = sweep.run_rules([_rule("class-assigned"), _rule("conditions-total")], everything)
    failures = {name: violations[:5] for name, violations in found.items() if violations}
    assert not failures, "\n".join(f"{name}: {len(found[name])} bodies, first {violations}" for name, violations in failures.items())


@pytest.mark.slow
def test_the_stored_class_is_the_conditions_class_or_the_anomaly_is_known(everything, expectations):
    """The generator still draws the class letter (GEN.646). While `class-from-conditions` is tagged the mismatches are reported;
    without the tag any mismatch fails. A tag on a rule that now holds everywhere is listed as stale in the report."""
    found = sweep.run_rules([_rule("class-from-conditions")], everything)["class-from-conditions"]
    tag = expectations["anomalies"].get("class-from-conditions")
    lines = [f"class-from-conditions: {len(found)} of {len(everything)} bodies; known anomaly {tag.split(':')[0]}" if tag
             else f"class-from-conditions: {len(found)} of {len(everything)} bodies"]
    lines += [f"  {label}: {violation}" for label, violation in found[:20]]
    if tag and not found:
        lines.append("STALE TAG: the rule holds everywhere; remove the anomaly from tests/data/class_sweep_expectations.json")
    census.write_report("\n".join(lines) + "\n", name=REPORT_NAME)
    assert tag or not found, f"{len(found)} bodies whose stored class is not the conditions' class, first {found[:5]}"


@pytest.mark.slow
def test_each_class_rule_reports_a_body_with_exactly_that_rule_broken(everything):
    assert sweep.prove_rules_can_fail(rules.RULES, everything) == {rule.name: True for rule in rules.RULES}


@pytest.mark.slow
def test_a_request_for_a_habitable_world_is_answered_or_the_anomaly_is_known(expectations):
    statuses = rules.habitable_request_statuses()
    share = rules.request_met_share(statuses)
    tag = expectations["anomalies"].get("habitable-request")
    minimum = expectations["habitable_request_minimum_met_share"]
    lines = [f"+habitable_world: {share:.0%} of {len(statuses)} systems answered aim or fallback ({statuses.count('unmet')} unmet)"
             + (f"; known anomaly {tag.split(':')[0]}" if tag else "")]
    if tag and share >= minimum:
        lines.append("STALE TAG: every request is answered; remove the anomaly from tests/data/class_sweep_expectations.json")
    census.write_report("\n".join(lines) + "\n", name="habitable-request-latest.md")
    assert tag or share >= minimum, f"only {share:.0%} of +habitable_world requests were answered: {statuses}"


@pytest.mark.slow
def test_the_same_seed_gives_the_same_classes_twice():
    assert rules.systems_with_a_class_from_conditions()
