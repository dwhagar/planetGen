# tests/test_score_filters.py

"""
UX.214: the search filters and sorts on the PHI-6 colour, the six domains, the four life scores and the plain life questions,
checked against the scores a generated sample galaxy (`tests/habitability_site.py`) stored: every filter must return exactly
the bodies a plain reading of the stored columns picks, and every sort must put them in order.
"""

import pytest

from planetgen.db import query, store
from planetgen.physics import habitability, score_filters
from planetgen.physics.habitability_world import DOMAINS, TIER_NAMES

from tests.habitability_site import habitable_galaxy  # noqa: F401 -- fixture

LIFE = ("microbial", "flora", "fauna", "sentient")


def _bodies(config, table):
    conn = store.get_connection(config)
    try:
        return [dict(row) for row in conn.execute(f"SELECT * FROM {table} WHERE phi4 IS NOT NULL").fetchall()]
    finally:
        conn.close()


def _found(config, panel, facet, value, sort=None):
    conn = query.open_readonly(config)
    try:
        tags = {name: set() for name in query.SEARCH_TAG_FACETS}
        tags[facet] = {value}
        result = query.search(conn, {}, tags, sorts={panel[:-1]: sort} if sort else None, limit=300)
        return result["results"][panel]["rows"] if result["results"][panel] else []
    finally:
        conn.close()


def _expected(row, facet, value):
    life = {name: row[f"life_{name}"] for name in LIFE}
    if facet == "habitability":
        return habitability.tier(row["phi4"]).lower() == value
    if facet == "domain":
        domain, colour = value.split(":")
        return TIER_NAMES[row[f"tier_{domain}"]].lower() == colour
    if facet == "life_level":
        name, least = value.split(":")
        return life[name] >= int(least)
    return {"alive": any(life.values()), "lifeless": not any(life.values()),
            "food": any(life[name] for name in ("flora", "fauna", "sentient")), "civilization": life["sentient"] > 0}[value]


@pytest.mark.parametrize("table,panel,prefix", (("planets", "planets", ""), ("moons", "moons", "moon_")))
def test_every_filter_returns_the_bodies_the_stored_scores_pick(habitable_galaxy, table, panel, prefix):
    bodies = _bodies(habitable_galaxy.config, table)
    assert len(bodies) > 3
    for facet, (_title, ticks) in score_filters.FACETS.items():
        for value, _label, _sql in ticks:
            want = {row["name"] for row in bodies if _expected(row, facet, value)}
            got = {row["name"] for row in _found(habitable_galaxy.config, panel, prefix + facet, value)}
            assert got == want, (facet, value, got ^ want)


def test_a_stray_value_matches_nothing(habitable_galaxy):
    assert _found(habitable_galaxy.config, "planets", "habitability", "purple") == []


@pytest.mark.parametrize("key,column", sorted(score_filters.SORT_COLUMNS.items()))
def test_a_sort_puts_the_worst_first_and_with_a_minus_the_best_first(habitable_galaxy, key, column):
    scores = lambda rows: [next(b[column] for b in _bodies(habitable_galaxy.config, "planets") if b["name"] == row["name"])
                           for row in rows]
    conn = query.open_readonly(habitable_galaxy.config)
    try:
        tags = {name: set() for name in query.SEARCH_TAG_FACETS}
        tags["living"] = {"alive", "lifeless"}          # every scored body, so each sort has something to order
        worst = query.search(conn, {}, tags, sorts={"planet": key}, limit=300)["results"]["planets"]["rows"]
        best = query.search(conn, {}, tags, sorts={"planet": f"-{key}"}, limit=300)["results"]["planets"]["rows"]
    finally:
        conn.close()
    assert worst and {r["name"] for r in worst} == {r["name"] for r in best}
    worst_scores, best_scores = scores(worst), scores(best)
    if key in score_filters.HIGH_FIRST:
        assert worst_scores == sorted(worst_scores) and best_scores == sorted(best_scores, reverse=True)
    else:       # a domain colour is stored 0 Blue to 3 Red, so worst (Red) first runs the number down
        assert worst_scores == sorted(worst_scores, reverse=True) and best_scores == sorted(best_scores)


def test_the_search_answers_the_options_each_filter_offers(habitable_galaxy):
    conn = query.open_readonly(habitable_galaxy.config)
    try:
        facets = query.search(conn, {}, {name: set() for name in query.SEARCH_TAG_FACETS})["facets"]
    finally:
        conn.close()
    for name in score_filters.FACETS:
        for prefix in ("", "moon_"):
            assert all(option["count"] > 0 for option in facets[prefix + name])
    assert {option["value"] for option in facets["habitability"]} <= score_filters.values("habitability")
    assert facets["habitability"], "the sample galaxy has no scored planet"


def test_the_chips_are_named_in_words(habitable_galaxy):
    conn = query.open_readonly(habitable_galaxy.config)
    try:
        answer = query.search(conn, {}, {name: set() for name in query.SEARCH_TAG_FACETS})
    finally:
        conn.close()
    labels = answer["facet_labels"]
    for key, label in (("living:alive", "Has life"), ("living:civilization", "Has a civilization")):
        if key in labels:
            assert labels[key] == label
    assert all(not label.startswith(("blue", "green")) for key, label in labels.items() if key.startswith("habitability:"))
    assert any(key.startswith("habitability:") and "PHI-6" in label for key, label in labels.items())
