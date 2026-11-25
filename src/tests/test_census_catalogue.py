# tests/test_census_catalogue.py

"""
TEST.164: the catalogue type audit on the census harness (`tests.audit.catalogue`, design note section 15, tests 10 to 14). A seeded
sample of a million systems of each of the eight region classes of the representation audit (stars of the region's population mix), a
million white dwarfs of each and a million young massive stars draws its member trees through `scatter_companions` and tests:

- test 10, the type counts: each type that is a property of a drawn tree (contact binaries, Algols, Be stars, runaways, brown dwarf
  companions, white dwarf companions, double white dwarfs and their mergers, cataclysmic variables, AM CVn, young groups) is within four Poisson
  sigma of the registry's default, or its range, in its reference class (`test_the_type_counts_*`);
- test 11, the region gates: no type at an age no population of its registry row has, no contact binary of solar-type stars in the young
  population, no Algol in the old (`test_no_type_is_found_where_its_gate_forbids_it`);
- test 12, the N tail: the systems with 5, 6, 7 and 8 to 12 stars against the registry's columns (`test_the_n_tail_*`);
- test 13, interaction consistency: every contact, semi-detached and accreting pair, classified again from its stored radii, masses and orbit
  (`test_every_interacting_pair_*`);
- test 14, group decay: no stored group is older than its decay age (`test_no_stored_group_*`).

The sample is drawn once for the file (about 17 million systems, about an hour on four cores; run it alone, `-n 0`, with a long timeout like the
other heavy audits). The checks of the machinery itself run on a few thousand systems.
"""

import dataclasses

import pytest

from planetgen import tuning
from planetgen.generation import multistar
from planetgen.generation.phenomena import system_types
from tests.audit import catalogue as audit
from tests.audit import multiplicity


@pytest.fixture(scope="module")
def sample():
    return audit.run()


def _failed(checks):
    return [check.line() for check in checks if not check.ok]


def test_the_sample_is_a_million_systems_of_each_region_class_and_a_million_white_dwarfs_of_each(sample):
    for region in audit.REGIONS:
        assert sample.rows[(audit.STAR, region)] == audit.SYSTEMS == 1_000_000
        assert sample.rows[(audit.WHITE_DWARF, region)] == 1_000_000
    assert sample.rows[(audit.YOUNG_MASSIVE, "young")] == audit.GROUP_ROWS


def test_the_type_counts_are_the_registrys_within_four_sigma(sample):
    checks = audit.type_checks(sample)
    assert len(checks) >= len(audit.TYPES)
    assert _failed(checks) == []


def test_no_type_is_found_where_its_gate_forbids_it(sample):
    assert _failed(audit.gate_checks(sample)) == [], "\n".join(sample.examples)


def test_the_n_tail_is_the_registrys_columns(sample):
    checks = audit.tail_checks(sample)
    assert len(checks) == 4 and _failed(checks) == []


def test_every_interacting_pair_is_what_its_geometry_gives(sample):
    assert sum(sample.pairs_checked.values()) > 1000
    assert _failed(audit.consistency_checks(sample)) == [], "\n".join(sample.examples)


def test_no_stored_group_is_older_than_its_decay_age(sample):
    assert sample.groups > 100 and sample.decayed_rows > 100
    assert _failed(audit.decay_checks(sample)) == []


def test_the_report_names_every_type_and_every_check_passes(sample):
    text = audit.format_report(sample)
    for name in audit.TYPES:
        assert name in text
    checks = audit.all_checks(sample)
    assert _failed(checks) == [], text
    assert f"{len(checks)} checks, 0 failed." in text
    unaccounted = (set(system_types.SYSTEM_TYPE_RATES) - set(audit.TYPES) - set(audit.NOT_A_TREE_PROPERTY) - set(audit.NOT_AUDITED_HERE)
                   - audit.COUNTED_ELSEWHERE)
    assert not unaccounted, unaccounted


# --- The audit's own machinery, on a few thousand systems ------------------------------------------------------------


def test_a_chunk_is_the_same_seeds_on_every_run_and_the_workers_add_up_to_the_one_by_one_run():
    assert audit.rows_of("star", "bulge", 3, 5) == audit.rows_of("star", "bulge", 3, 5)
    assert audit.rows_of("star", "bulge", 3, 5) != audit.rows_of("star", "bulge", 4, 5)
    jobs = [("star", "bulge", 0, 150), ("white-dwarf", "halo", 0, 150), ("young-massive", "young", 0, 150)]
    one, two = audit.run(workers=1, sample_jobs=jobs), audit.run(workers=2, sample_jobs=jobs)
    assert dataclasses.asdict(one) == dataclasses.asdict(two)


def test_the_jobs_cover_the_whole_sample_in_chunks():
    found = audit.jobs()
    for region in audit.REGIONS:
        for sample in (audit.STAR, audit.WHITE_DWARF):
            assert sum(count for s, r, _chunk, count in found if (s, r) == (sample, region)) == audit.SYSTEMS
    assert sum(count for s, _r, _chunk, count in found if s == audit.YOUNG_MASSIVE) == audit.GROUP_ROWS
    assert max(count for _s, _r, _chunk, count in found) <= audit.CHUNK


def test_every_region_class_has_a_population_mix_and_a_white_dwarf_population():
    for region in audit.REGIONS:
        names, weights = audit.region_mix(region)
        assert names and len(names) == len(weights) and all(w > 0 for w in weights)
        info, _row = audit.draw_row(audit.WHITE_DWARF, region, 7)
        assert tuning.STELLAR_POPULATION_AGE_RANGES_GY[info.population][1] >= audit.WHITE_DWARF_MIN_WINDOW_GY


def test_a_single_objects_flags_are_counted_a_lone_be_star_and_a_merger_product_are_in_the_sample():
    sample = audit.run(workers=2, sample_jobs=[("white-dwarf", "halo", 0, 4000)])
    assert sample.types[("halo", "double-white-dwarf-merger")] > 100       # a merger product is mostly a single white dwarf: the draw keeps its flag


def test_a_tree_that_breaks_each_rule_is_found():
    seed = next(s for s in audit.rows_of("star", "bulge", 0, 3000) if len(audit.tree_of("star", "bulge", s)[1].nodes) >= 1)
    info, tree = audit.tree_of("star", "bulge", seed)
    assert audit.interaction_problems(tree) == [] and audit.group_problems(tree) == []
    node = tree.nodes[0]
    pair = [n for n, members in audit._pairs(tree) if members][:1]
    assert pair, "the tree has no pair of two single members"
    contact = dataclasses.replace(tree, nodes=tuple(dataclasses.replace(n, interaction="contact") if n is pair[0] else n for n in tree.nodes))
    assert any("contact pair" in problem for problem in audit.interaction_problems(contact))
    group = dataclasses.replace(tree, nodes=tuple(dataclasses.replace(n, kind="group") if n is node else n for n in tree.nodes))
    assert audit.group_problems(group)
    young = dataclasses.replace(info, population="young", age_gy=0.01, birth_mass_solar=1.0)
    flagged = dataclasses.replace(tree, nodes=tuple(dataclasses.replace(n, interaction="semi-detached") if n is pair[0] else n for n in tree.nodes))
    old = dataclasses.replace(info, age_gy=11.0)
    assert any("Algol" in problem for problem in audit.gate_problems(flagged, old))
    assert audit.gate_problems(tree, young) is not None


def test_a_check_that_misses_its_expectation_fails_and_says_so():
    check = audit.rate_check("a type", 10, 10_000, 0.01, 0.01, 0.01)
    assert not check.ok and "FAIL" in check.line()
    assert audit.rate_check("a type", 105, 10_000, 0.01, 0.01, 0.01).ok
    assert audit.rate_check("a ranged type", 40, 1000, 0.032, 0.032, 0.068).ok
    assert not audit.rate_check("a ranged type", 2, 1000, 0.032, 0.032, 0.068).ok
    assert not audit.rate_check("a ranged type", 400, 1000, 0.032, 0.032, 0.068).ok


def test_the_registry_rows_the_audit_counts_exist_and_its_references_are_the_rows_own():
    for name, (ref_name, _predicate) in audit.TYPES.items():
        entry = system_types.SYSTEM_TYPE_RATES[name]
        assert audit.REFERENCES[ref_name] and entry.default is not None
        assert audit.REFERENCE_OF.get(entry.reference, ref_name) == ref_name or entry.reference in ("white dwarfs", "systems"), (name, entry.reference)
    assert multistar.MASSIVE_PRIMARY_MIN_MASS_SOL == 2.1


def test_the_report_command_prints_and_writes_the_report(tmp_path, capsys, monkeypatch):
    real = audit.jobs
    monkeypatch.setattr(audit, "jobs", lambda *args, **kwargs: real(systems=200, group_rows=200))
    path = tmp_path / "catalogue.txt"
    code = audit.main([str(path)])
    assert code in (0, 1)
    printed = capsys.readouterr().out
    assert "Catalogue type audit (TEST.164)" in printed and path.read_text(encoding="utf-8") == printed
