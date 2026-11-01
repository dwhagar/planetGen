# tests/test_host_tables.py

"""GEN.552: the host classes and their gated tables. The tables only move what the code did into rows, so every draw is
checked against the formula it replaced (`random() < odds`, `choices(range(1, 8), weights)`) over many seeds."""

from dataclasses import replace

import pytest

from planetgen.generation import possibility_table as pt
from planetgen.generation import system_pass
from planetgen.generation.phenomena import hosts, registry
from planetgen.util import draw

ORIGINAL_MILLISECOND_PULSAR_PLANETS = 0.007
ORIGINAL_ROGUE_MOON = 0.2


def host_of(kind):
    """The `HostClass` of the registry kind `kind`, or `None` for a kind that hosts nothing (a nebula, a comet)."""
    name = registry.BY_KIND[kind].host
    return None if name is None else hosts.BY_NAME[name]


def reweighted(table, row_name, weight):
    """A copy of `table` with the row `row_name` given `weight`; the other tests use it to move one row of a host table."""
    rows = [replace(row, weight=weight) if row.name == row_name else row for row in table.rows]
    return pt.Table(table.name, table.version, rows, table.components, table.description)


def test_every_kind_names_a_host_class_or_none_and_every_class_has_its_tables():
    for row in registry.KINDS:
        assert row.host is None or row.host in hosts.BY_NAME, row.kind
    assert {row.host for row in registry.KINDS} >= {"star-like", "remnant", "planet", "brown-dwarf"}
    for host in hosts.HOSTS:
        for name in (host.companions, host.satellites):
            assert name is None or name in pt.REGISTRY, (host.name, name)
    assert hosts.BY_NAME["moon"].companions is None and hosts.BY_NAME["moon"].satellites == "moon-satellites"
    assert host_of("nebula") is None and host_of("neutron-star").name == "remnant"


def test_a_moon_and_a_small_body_have_host_rows_gated_by_may_host():
    """Boss's cards S5 and S6: moons are no longer the lowest level. A small body's gate is closed and it has no table; a
    moon's table has only the empty row until GEN.579's tidal-survival row, so neither holds anything today."""
    small = hosts.BY_NAME["small-body"]
    assert not small.may_host and small.satellites is None
    assert host_of("comet") is small and host_of("asteroid-field") is small
    moon = hosts.BY_NAME["moon"]
    assert moon.may_host
    assert [row.name for row in hosts.REGISTRY_TABLES["moon-satellites"].rows] == [hosts.NO_SATELLITES]
    assert not hosts.holds("moon-satellites", "moon") and hosts.chance("moon-satellites", "moon") == 0.0
    assert all(host.satellites is None or host.may_host for host in hosts.HOSTS)


def test_the_tables_are_registered_and_each_row_says_where_its_weight_comes_from():
    for name in ("star-count", "remnant-satellites", "planet-satellites", "moon-satellites"):
        table = pt.REGISTRY.get(name)
        assert table is hosts.REGISTRY_TABLES[name]
        for row in table.rows:
            assert row.source and row.state in pt.STATES, (name, row.name)
    assert "stand-in" in pt.REGISTRY.get("planet-satellites").state_counts()


def test_the_star_count_table_holds_the_class_weights_and_scales_the_five_or_more_columns():
    for letter, weights in hosts.STAR_COUNT_PROBABILITIES.items():
        assert sum(weights) == pytest.approx(1.0, abs=1e-4), letter
        assert hosts.star_count_weights(letter) == pytest.approx(weights)
        scaled = hosts.star_count_weights(letter, 0.5)
        assert scaled[:4] == pytest.approx(weights[:4]) and scaled[4:] == pytest.approx([w * 0.5 for w in weights[4:]])


@pytest.mark.parametrize("letter, prevalence", [("O", 1.0), ("G", 1.0), ("B-mid", 1.0), ("M", 1.0), ("O", 0.0),
                                                ("O", 2.5), ("K", 0.3)])
def test_the_star_count_draw_is_the_original_formula_on_the_same_stream(letter, prevalence):
    weights = tuple(w * (prevalence if count >= 5 else 1.0)
                    for count, w in enumerate(hosts.STAR_COUNT_PROBABILITIES[letter], start=1))
    for seed in range(1500):
        original = draw.Stream(seed).choices(range(1, 8), weights)[0]
        rng = draw.Stream(seed)
        assert hosts.draw_star_count(letter, prevalence, rng) == original, (seed, letter, prevalence)
        # The same single draw leaves the stream where the original left it.
        original_stream = draw.Stream(seed)
        original_stream.choices(range(1, 8), weights)
        assert rng.random() == original_stream.random()


def test_a_millisecond_pulsar_rolls_planets_exactly_when_the_original_odds_would():
    key = "neutron-star:millisecond"
    for seed in range(30_000):
        original = draw.Stream(f"{seed}:{system_pass.STREAM}").random() < ORIGINAL_MILLISECOND_PULSAR_PLANETS
        assert system_pass.hit(key, seed) == original, seed


def test_the_other_remnants_and_kinds_hold_nothing():
    for key in ("neutron-star:magnetar", "neutron-star:young", "black-hole:stellar", "black-hole:intermediate",
                "black-hole:supermassive", "comet", "asteroid-field", "brown-dwarf", "rogue-planet"):
        assert not hosts.holds("remnant-satellites", key), key
        assert hosts.chance("remnant-satellites", key) == 0.0


def test_a_rogue_planet_keeps_a_moon_exactly_when_the_original_odds_would():
    for seed in range(30_000):
        original = draw.Stream(seed).random() < ORIGINAL_ROGUE_MOON
        rng = draw.Stream(seed)
        assert hosts.holds_something("planet-satellites", "rogue-planet", rng) == original, seed
        original_stream = draw.Stream(seed)
        original_stream.random()
        assert rng.random() == original_stream.random()


def test_the_no_satellites_row_is_the_rest_of_the_chance():
    assert hosts.chance("remnant-satellites", "neutron-star:millisecond") == pytest.approx(0.007)
    probabilities = dict((row.name, p) for row, p in hosts.REMNANT_SATELLITES.probabilities(
        pt.Subject("neutron-star:millisecond")))
    assert probabilities == {"planets": pytest.approx(0.007), hosts.NO_SATELLITES: pytest.approx(0.993)}
    assert dict((row.name, p) for row, p in hosts.REMNANT_SATELLITES.probabilities(pt.Subject("black-hole:stellar"))) \
        == {hosts.NO_SATELLITES: 1.0}


def test_a_reweighted_table_changes_only_the_named_row():
    table = reweighted(hosts.REMNANT_SATELLITES, "planets", 0.5)
    assert table.row("planets").weight == 0.5
    assert table.row(hosts.NO_SATELLITES) == hosts.REMNANT_SATELLITES.row(hosts.NO_SATELLITES)
    assert table.digest() != hosts.REMNANT_SATELLITES.digest()
