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

from tests.scored_site import scored_galaxy  # noqa: F401 -- fixture

LIFE = ("microbial", "flora", "fauna", "sentient")


def _bodies(config, table):
    conn = store.get_connection(config)
    try:
        return [dict(row) for row in conn.execute(f"SELECT * FROM {table} WHERE phi4 IS NOT NULL").fetchall()]
    finally:
        conn.close()


def _class_scores(config, table):
    """`{body name: {class id: score}}` from the body's life class score table (MAP.345)."""
    kind = table[:-1]
    scores, owner = score_filters.SCORE_TABLES[kind]
    conn = store.get_connection(config)
    try:
        found = {}
        for row in conn.execute(f"SELECT b.name, s.life_class, s.score FROM {scores} s JOIN {table} b ON b.id = s.{owner}").fetchall():
            found.setdefault(row["name"], {})[row["life_class"]] = row["score"]
        return found
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


def _expected(row, facet, value, class_scores=None):
    if facet == "supports":
        class_id, least = value.split(":")
        return (class_scores or {}).get(row["name"], {}).get(class_id, 0) >= int(least)
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
def test_every_filter_returns_the_bodies_the_stored_scores_pick(scored_galaxy, table, panel, prefix):
    bodies = _bodies(scored_galaxy.config, table)
    class_scores = _class_scores(scored_galaxy.config, table)
    assert len(bodies) > 3
    for facet, (_title, ticks) in score_filters.FACETS.items():
        for value, _label, _sql in ticks:
            want = {row["name"] for row in bodies if _expected(row, facet, value, class_scores)}
            got = {row["name"] for row in _found(scored_galaxy.config, panel, prefix + facet, value)}
            assert got == want, (facet, value, got ^ want)


def test_a_stray_value_matches_nothing(scored_galaxy):
    assert _found(scored_galaxy.config, "planets", "habitability", "purple") == []


@pytest.mark.parametrize("key,column", sorted(score_filters.SORT_COLUMNS.items()))
def test_a_sort_puts_the_worst_first_and_with_a_minus_the_best_first(scored_galaxy, key, column):
    scores = lambda rows: [next(b[column] for b in _bodies(scored_galaxy.config, "planets") if b["name"] == row["name"])
                           for row in rows]
    conn = query.open_readonly(scored_galaxy.config)
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


def test_the_search_answers_the_options_each_filter_offers(scored_galaxy):
    conn = query.open_readonly(scored_galaxy.config)
    try:
        facets = query.search(conn, {}, {name: set() for name in query.SEARCH_TAG_FACETS})["facets"]
    finally:
        conn.close()
    for name in score_filters.FACETS:
        for prefix in ("", "moon_"):
            assert all(option["count"] > 0 for option in facets[prefix + name])
    assert {option["value"] for option in facets["habitability"]} <= score_filters.values("habitability")
    assert facets["habitability"], "the sample galaxy has no scored planet"


def test_the_chips_are_named_in_words(scored_galaxy):
    conn = query.open_readonly(scored_galaxy.config)
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


@pytest.mark.parametrize("table,prefix", (("planets", ""), ("moons", "moon_")))
def test_the_supports_options_count_the_bodies_at_each_class_and_threshold(scored_galaxy, table, prefix):
    class_scores = _class_scores(scored_galaxy.config, table)
    conn = query.open_readonly(scored_galaxy.config)
    try:
        facets = query.search(conn, {}, {name: set() for name in query.SEARCH_TAG_FACETS})["facets"]
    finally:
        conn.close()
    offered = {option["value"]: option["count"] for option in facets[prefix + "supports"]}
    for value, _label, _sql in score_filters.FACETS["supports"][1]:
        class_id, least = value.split(":")
        want = sum(1 for scores in class_scores.values() if scores.get(class_id, 0) >= int(least))
        assert offered.get(value, 0) == want, value
    assert (class_scores and offered) or table == "moons", "the sample galaxy stored no class scores, or none is offered"


def test_a_body_that_supports_a_class_is_found_by_the_class_and_not_by_another(scored_galaxy):
    class_scores = _class_scores(scored_galaxy.config, "planets")
    supported = {(name, class_id) for name, scores in class_scores.items() for class_id in scores}
    if not supported:
        pytest.skip("the sample galaxy stored no planet class scores")
    name, class_id = sorted(supported)[0]
    assert name in {row["name"] for row in _found(scored_galaxy.config, "planets", "supports", f"{class_id}:1")}
    others = [c.id for c in __import__("planetgen.physics.life_classes", fromlist=["CLASSES"]).CLASSES
              if (name, c.id) not in supported]
    assert others and name not in {row["name"] for row in _found(scored_galaxy.config, "planets", "supports", f"{others[0]}:1")}
