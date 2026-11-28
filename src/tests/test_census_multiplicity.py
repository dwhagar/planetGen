# tests/test_census_multiplicity.py

"""
TEST.162: the multiplicity audit on the census harness (`tests.audit.multiplicity`). A seeded scatter of a million rows of each placed
class of the phenomenon scatter (a stellar and an intermediate-mass black hole, a neutron star, a planetary nebula's central star), a
million star rows over the nine spectral classes and a million white dwarf rows draw their member trees through the function the scatter
passes call, and every row and tree is tested:

- the share of rows with each companion state is the table's within four Poisson sigma, and the stars of a star row follow the
  `star-count` table (`test_the_companion_states_*`, `test_the_stars_of_a_star_row_*`);
- a class's rows and their same-class companions keep the class's object density; nobody is both a row and a member
  (`test_the_object_density_*`, `test_nobody_is_both_a_row_and_a_member`);
- every stored tree passes the stability, Kozai and Peters screens, the same seed draws the same members, and nothing in a tree is a
  planet, a moon or a belt (`test_every_stored_tree_*`, `test_the_same_seed_*`, `test_nothing_in_a_scatter_is_a_planet_*`);
- black hole pairs lie between the gravitational-wave stock and the Lamberts count (`test_black_hole_pairs_*`).

The database checks (a built system equals its stored members, the nightly update writes no member row) are in
`test_census_multiplicity_db.py`; the report of this file is also a command, `python -m tests.audit.multiplicity [report.txt]`, and
`python -m tests.audit.representation` prints it with the representation audit's.

The sample is drawn once for the file (about 1.4 million trees; run it alone, `-n 0`, with a long timeout like the other heavy audits).
"""

import dataclasses

import pytest

from planetgen.generation import multistar
from tests.audit import multiplicity as audit


@pytest.fixture(scope="module")
def sample():
    return audit.run()


def _failed(checks):
    return [check.line() for check in checks if not check.ok]


def test_the_sample_is_at_least_a_million_rows_of_each_class(sample):
    for cls in audit.PHENOMENON_CLASSES:
        assert sample[cls.name].rows == audit.ROWS == 1_000_000
    assert sum(sample[name].rows for name, _lo, _hi in audit.STAR_CLASSES) >= 1_000_000
    assert sample[audit.WHITE_DWARF].rows == 1_000_000
    assert set(sample) == set(audit.classes())


def test_the_companion_states_of_every_phenomenon_class_are_the_tables_within_four_sigma(sample):
    checks = audit.companion_state_checks(sample)
    assert len(checks) >= 4 * 4
    assert _failed(checks) == []


def test_the_stars_of_a_star_row_follow_the_star_count_table_class_by_class(sample):
    checks = audit.star_count_checks(sample)
    assert len(checks) == 3 * len(audit.STAR_CLASSES)
    assert _failed(checks) == []


def test_a_white_dwarf_row_has_the_companion_odds_of_the_white_dwarf_table(sample):
    tally = sample[audit.WHITE_DWARF]
    table = audit.state_table_of("white-dwarf")
    for state, odds in table.items():
        check = audit.share_check(f"white-dwarf rows with a {state} companion", tally.states[state], tally.rows, odds)
        assert check.ok, check.line()
    none = audit.share_check("white-dwarf rows with no companion", tally.states[None], tally.rows, 1.0 - sum(table.values()))
    assert none.ok, none.line()


def test_the_object_density_of_a_class_is_the_same_with_rows_and_members(sample):
    checks = [c for c in audit.bookkeeping_checks(sample) if "same-class" in c.name or "divisor" in c.name]
    assert checks and _failed(checks) == []


def test_nobody_is_both_a_row_and_a_member(sample):
    for name in audit.classes():
        assert sample[name].both_row_and_member == 0, name
        assert sample[name].duplicate_member_seeds == 0, name


def test_every_stored_tree_passes_the_stability_kozai_and_peters_screens(sample):
    checks = [c for c in audit.tree_checks(sample) if "failing the stability" in c.name]
    assert len(checks) == len(audit.classes())
    assert _failed(checks) == [], "\n".join(example for name in audit.classes() for example in sample[name].unstable_examples)
    assert sum(sample[name].trees for name in audit.classes()) > 1_000_000   # the screens were asked of over a million trees


def test_the_same_seed_draws_the_same_members(sample):
    for name in audit.classes():
        assert sample[name].determinism_checked > 0, name
        assert sample[name].determinism_broken == 0, name


def test_nothing_in_a_scatter_is_a_planet_a_moon_or_a_belt(sample):
    for name in audit.classes():
        assert not sample[name].not_a_member_kind, (name, dict(sample[name].not_a_member_kind))
    assert set(multistar.MEMBER_KINDS) == {"main-sequence", "giant", "white-dwarf", "neutron-star", "black-hole", "brown-dwarf",
                                           "pre-main-sequence", "protostar", "stripped-helium", "thorne-zytkow", "wolf-rayet"}


def test_black_hole_pairs_lie_between_the_gravitational_wave_stock_and_the_lamberts_count(sample):
    tally = sample["black-hole"]
    assert tally.rows == 1_000_000 and tally.black_hole_pairs > 0
    share = tally.black_hole_pairs / tally.rows
    low, high = audit.BLACK_HOLE_PAIR_SHARE
    assert low <= share <= high, share
    assert _failed(audit.black_hole_pair_checks(sample)) == []


def test_the_report_names_every_class_and_every_check_passes(sample):
    text = audit.format_report(sample)
    for name in audit.classes():
        assert name in text
    checks = audit.all_checks(sample)
    assert _failed(checks) == [], text
    assert f"{len(checks)} checks, 0 failed." in text


# --- The audit's own machinery, on a few thousand rows --------------------------------------------------------------


SMALL = [(name, 0, 2000) for name in audit.classes()]


def test_a_chunk_is_the_same_seeds_on_every_run_and_the_workers_add_up_to_the_one_by_one_run():
    assert audit.rows_of("black-hole", 3, 5) == audit.rows_of("black-hole", 3, 5)
    assert audit.rows_of("black-hole", 3, 5) != audit.rows_of("black-hole", 4, 5)
    jobs = [("black-hole", 0, 3000), ("G", 0, 300), ("white-dwarf", 1, 300)]
    one = audit.run(workers=1, sample_jobs=jobs)
    two = audit.run(workers=2, sample_jobs=jobs)
    assert {name: dataclasses.asdict(tally) for name, tally in one.items()} == {name: dataclasses.asdict(tally) for name, tally in two.items()}


def test_the_jobs_cover_the_whole_sample_in_chunks():
    found = audit.jobs()
    for name in audit.classes():
        assert sum(count for job_name, _chunk, count in found if job_name == name) == audit.rows_in(name)
    assert max(count for _name, _chunk, count in found) <= audit.CHUNK


def test_an_unstable_tree_is_found_by_the_stability_screen():
    seed = next(s for s in audit.rows_of("O", 0, 2000) if (audit.tree_of("O", s)[1] is not None and audit.tree_of("O", s)[1].member_count >= 3))
    tree = audit.tree_of("O", seed)[1]
    assert audit.stability_problems(tree) == []
    nodes = list(tree.nodes)
    inner = next(node for node in nodes if node.parent_node_no is not None)
    squeezed = [dataclasses.replace(n, separation_au=n.separation_au * 1.0e3, separation_km=n.separation_km * 1.0e3) if n is inner else n
                for n in nodes]
    assert audit.stability_problems(dataclasses.replace(tree, nodes=tuple(squeezed)))


def test_a_check_that_misses_its_expectation_fails_and_says_so():
    check = audit.share_check("a share", 1000, 10_000, 0.2)
    assert not check.ok and "FAIL" in check.line()
    assert audit.share_check("a share", 2010, 10_000, 0.2).ok
    assert audit.Check("one side", 5.0, 3.0, 1.0, 4.0, "at most").ok
    assert not audit.Check("one side", 9.0, 3.0, 1.0, 4.0, "at most").ok
    assert audit.Check("one side", 1.0, 3.0, 1.0, 4.0, "at least").ok
    assert not audit.Check("one side", -2.0, 3.0, 1.0, 4.0, "at least").ok


def test_the_report_command_prints_and_writes_the_report(tmp_path, capsys, monkeypatch):
    jobs = [(name, 0, 500) for name in audit.classes()]
    monkeypatch.setattr(audit, "jobs", lambda: jobs)
    path = tmp_path / "report.txt"
    code = audit.main([str(path)])
    assert code in (0, 1)
    printed = capsys.readouterr().out
    assert "Multiplicity audit (TEST.162)" in printed and path.read_text(encoding="utf-8") == printed
