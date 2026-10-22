# tests/test_census_representation.py

"""
TEST.160: the representation audit. Seeded sectors are filled in each region class of the galaxy (the bulge, which holds
the bar; a spiral arm and the space between arms; the thick disk; the halo; a Magellanic-type satellite; a globular
cluster; the intergalactic floor), counted by stage and kind on the census harness (`tests.audit.census`), and compared
with what the profile's rates, stars and factors expect: every kind within four Poisson sigma, a kind expected 20 times
or more appears at least once, and none appears where its gate forbids it. "All types represented" means audited
frequency, never forced (Boss, 2026-10-15). The rates are dialed up in the fills (`tests.audit.representation`).

What each part of the TODO item's text is covered by:

- counts per region and kind, the 20-times rule and the forbidden kinds: this file;
- every hosted nebula, planetary nebula, core-collapse remnant, D or E region and Wolf-Rayet nebula holds its required
  contents: `test_census_nebula_wholesphere.py` (TEST.184) and `test_census_component_wholesphere.py` (TEST.185);
- the share of systems with planets or a companion matches the system table: `test_gen_fill_stages.py`, `test_system_pass.py`.

The report is also a command: `python -m tests.audit.representation [report.txt]`.
"""

import pytest

from tests.audit import census as audit
from tests.audit import representation as rep


@pytest.fixture(scope="module")
def audit_run():
    patch = pytest.MonkeyPatch()
    try:
        yield rep.run(patch)
    finally:
        patch.undo()


def test_the_audit_covers_every_region_class_of_the_galaxy(audit_run):
    census, _fill, names = audit_run
    assert set(names) == {"bulge", "thin-disk-arm", "thin-disk-interarm", "thick-disk", "halo", "satellite",
                          "globular-cluster", "intergalactic"}
    for name in names:
        assert census[name].sectors == rep.SECTORS
        assert rep.dominant_population(name) == rep.REGION_POINTS[name][1]


def test_every_kind_is_within_four_sigma_of_its_expectation_in_every_region(audit_run):
    census, fill, names = audit_run
    expected = rep.expectations(fill, names)
    findings = audit.compare(census, expected)
    assert findings == [], audit.format_report(census, expected, findings)


def test_the_kinds_that_need_a_spot_of_their_own_never_exceed_their_rate(audit_run):
    census, fill, names = audit_run
    assert rep.limited_findings(census, fill, names) == []


def test_a_kind_expected_twenty_times_appears_and_the_audit_has_kinds_that_do(audit_run):
    census, fill, names = audit_run
    expected = rep.expectations(fill, names)
    appearing = [(name, kind) for name in names for kind, mean in expected[name].items()
                 if audit.expected_total(mean, census[name].sectors) >= audit.APPEAR_AT]
    assert len(appearing) >= 6
    for name, kind in appearing:
        assert census[name].observed(kind) > 0, (name, kind)


def test_no_kind_appears_where_the_places_gate_forbids_it(audit_run):
    census, _fill, _names = audit_run
    for name in ("halo", "globular-cluster", "intergalactic", "thick-disk"):
        for kind in ("molecular-cloud", "supernova-remnant", "nebula"):
            assert census[name].observed(kind) == 0, (name, kind)
    for name in ("thin-disk-arm", "satellite"):
        assert census[name].observed("molecular-cloud") > 0 and census[name].observed("supernova-remnant") > 0


def test_the_arm_shows_more_gas_than_the_space_between_arms(audit_run):
    census, _fill, _names = audit_run
    arm, between = rep.factors_at("thin-disk-arm"), rep.factors_at("thin-disk-interarm")
    assert arm["molecular-cloud"] > 10.0 * between["molecular-cloud"]
    assert census["thin-disk-arm"].observed("molecular-cloud") > census["thin-disk-interarm"].observed("molecular-cloud")


def test_the_report_names_every_region_and_ends_clean(audit_run):
    census, fill, names = audit_run
    expected = rep.expectations(fill, names)
    text = audit.format_report(census, expected, audit.compare(census, expected))
    for name in names:
        assert f"{name}: {rep.SECTORS} sectors" in text
    assert text.rstrip().endswith("Every expected count is within bounds")


def test_the_report_command_prints_and_writes_the_report(tmp_path, capsys):
    path = tmp_path / "representation.txt"
    assert rep.main([str(path), "--sectors", "6"]) in (0, 1)
    printed = capsys.readouterr().out
    written = path.read_text(encoding="utf-8")
    assert written.startswith("Census audit") and written in printed       # the fills' own log lines come before it
    assert "Placement-limited kinds" in written
