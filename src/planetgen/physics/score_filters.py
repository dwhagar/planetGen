# planetgen/physics/score_filters.py

"""
The habitability and life score filters of the search page and the API (UX.214): the colour of a body's PHI-6 score, the
colour of each of its six domains (Suitability among them), how much microbial, plant, animal and sentient life it has,
and the plain questions "is there life", "is there food" and "is there a civilization".

Each filter is a tick, as the survival equipment ones are (`equipment.TICKS`): `(value, label, SQL condition)`, the
condition on a `planets` or `moons` row with the table alias in `{a}`. A body matches a facet when any of its ticked
values holds. A body that was not scored (a gas giant, or one saved before the scores existed) matches none of them.
The overall colour reads the stored `phi4` through the same cut-offs as `habitability.tier`; the domains read the stored
`tier_<domain>` numbers (0 Blue to 3 Red); the life scores are whole numbers 0 to 9 (`life_scores`).
"""

from planetgen.physics import life_classes
from planetgen.physics.habitability_world import DOMAINS, TIER_NAMES

HABITABILITY_TICKS = (
    ("blue", "Blue (PHI-6 1 or more)", "{a}.phi4 >= 1.0"),
    ("green", "Green (PHI-6 0.75 to under 1)", "({a}.phi4 >= 0.75 AND {a}.phi4 < 1.0)"),
    ("yellow", "Yellow (PHI-6 0.4 to under 0.75)", "({a}.phi4 >= 0.4 AND {a}.phi4 < 0.75)"),
    ("red", "Red (PHI-6 under 0.4)", "{a}.phi4 < 0.4"),
)
"""tuple: The overall PHI-6 colour ticks."""

DOMAIN_TICKS = tuple(
    (f"{domain}:{name.lower()}", f"{domain.capitalize()}: {name}", f"{{a}}.tier_{domain} = {number}")
    for domain in DOMAINS for number, name in enumerate(TIER_NAMES))
"""tuple: One tick per domain and colour, `<domain>:<colour>` (for example `suitability:green`)."""

LIFE_KINDS = (("microbial", "Microbial life"), ("flora", "Flora"), ("fauna", "Fauna"), ("sentient", "Sentient life"))
LIFE_LEVELS = ((1, "present"), (5, "rich"))
"""tuple: The life score thresholds offered: 1 or more of 9 (present) and 5 or more of 9 (rich)."""

LIFE_LEVEL_TICKS = tuple(
    (f"{kind}:{level}", f"{label} {level} or more of 9 ({word})", f"{{a}}.life_{kind} >= {level}")
    for kind, label in LIFE_KINDS for level, word in LIFE_LEVELS)
"""tuple: One tick per life score and threshold, `<kind>:<least score>` (for example `flora:5`)."""

LIVING_TICKS = (
    ("alive", "Has life", "({a}.life_microbial > 0 OR {a}.life_flora > 0 OR {a}.life_fauna > 0 OR {a}.life_sentient > 0)"),
    ("lifeless", "Has no life",
     "({a}.life_microbial = 0 AND {a}.life_flora = 0 AND {a}.life_fauna = 0 AND {a}.life_sentient = 0)"),
    ("food", "Has food (plants, animals or a civilization)",
     "({a}.life_flora > 0 OR {a}.life_fauna > 0 OR {a}.life_sentient > 0)"),
    ("civilization", "Has a civilization", "{a}.life_sentient > 0"),
)
"""tuple: The plain life questions. `food` follows `habitability.food_score`: flora, fauna or a civilization is food."""

SUPPORTS_LEVELS = ((1, "supported"), (5, "well supported"))
"""tuple: The class score thresholds offered: 1 or more of 9 (the body supports the class at all) and 5 or more of 9."""

SUPPORTS_SQL = ("EXISTS (SELECT 1 FROM {scores} s WHERE s.{owner} = {a}.id AND s.life_class = '%(class_id)s' AND s.score >= %(least)d)")
"""str: The test of one "supports Class N" tick on a `planets` or `moons` row with the alias in `{a}`: the body has a row for the class in
its life class score table (`{scores}`, keyed by `{owner}`) at the least score. Only classes that reach the line have a row
(`life_class_scores`, GEN.889), so the tick reads the index on the class and the score."""

SUPPORTS_TICKS = tuple(
    (f"{definition.id}:{least}", f"Supports Class {definition.id}: {definition.name}, score {least} or more of 9 ({word})",
     SUPPORTS_SQL % {"class_id": definition.id, "least": least})
    for definition in life_classes.CLASSES for least, word in SUPPORTS_LEVELS)
"""tuple: One tick per life class and threshold, `<class id>:<least score>` (for example `1a:5`), "this body supports Class N life" with
the score the class got from the body's conditions (docs/design/habitability/life-classes.md section 10, MAP.345). The class that
lives on a body is the `life` facet."""

FACETS = {
    "habitability": ("Habitability colour", HABITABILITY_TICKS),
    "domain": ("Habitability domain colour", DOMAIN_TICKS),
    "life_level": ("Life score", LIFE_LEVEL_TICKS),
    "living": ("Life", LIVING_TICKS),
    "supports": ("Supports a life class", SUPPORTS_TICKS),
}
"""dict: Each facet (the search tag facet name without a `moon_` prefix) to its title and ticks."""

SCORE_TABLES = {"planet": ("planet_life_class_scores", "planet_id"), "moon": ("moon_life_class_scores", "moon_id")}
"""dict: The kind of a body to the table of its life class scores and the column that names the body there."""

SORT_COLUMNS = {
    "phi4": "phi4", "suitability": "tier_suitability", "microbial": "life_microbial", "flora": "life_flora",
    "fauna": "life_fauna", "sentient": "life_sentient",
    **{domain: f"tier_{domain}" for domain in DOMAINS},
}
"""dict: What a result list can be sorted by (`sort=<key>`), to the stored column."""

HIGH_FIRST = frozenset({"phi4", "microbial", "flora", "fauna", "sentient"})
"""frozenset: The sort keys whose stored number is higher for a better body (PHI-6 and the life scores). A domain colour is
stored 0 for Blue to 3 for Red, so it runs the other way."""


def values(facet):
    """The tick values of a facet."""
    return frozenset(value for value, _label, _sql in FACETS[facet][1])


def _owner(kind):
    scores, owner = SCORE_TABLES[kind]
    return {"scores": scores, "owner": owner}


def condition(facet, alias, selected, kind="planet"):
    """The SQL condition (a parenthesised `OR` of the selected ticks of `facet` on the table `alias`; unknown values are
    dropped) or `None` when none of `selected` is a tick. `kind` (`planet` or `moon`) names the body's life class score table
    for the `supports` facet."""
    conditions = [sql.format(a=alias, **_owner(kind)) for value, _label, sql in FACETS[facet][1] if value in selected]
    return "(" + " OR ".join(conditions) + ")" if conditions else None


def counts_select(facet, alias, kind="planet"):
    """The select list that counts the bodies of table `alias` each tick of `facet` holds for (`t0`, `t1`, ...)."""
    return ", ".join(f"COALESCE(SUM({sql.format(a=alias, **_owner(kind))}), 0) AS t{i}"
                     for i, (_v, _l, sql) in enumerate(FACETS[facet][1]))


def supports_counts_sql(kind):
    """The query that counts, for every class in the score table of `kind`, the bodies at each threshold of `SUPPORTS_LEVELS`
    (columns `life_class`, `t<least>`). One pass over the index on the class and the score, not a test for each body and tick."""
    scores, _owner_column = SCORE_TABLES[kind]
    sums = ", ".join(f"COALESCE(SUM(score >= {least}), 0) AS t{least}" for least, _word in SUPPORTS_LEVELS)
    return f"SELECT life_class, {sums} FROM {scores} GROUP BY life_class"
