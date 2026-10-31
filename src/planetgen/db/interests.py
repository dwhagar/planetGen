"""
The reasons a star system is worth a mark on the Galaxy Map (MAP.305, MAP.303): a registry of small SQL tests over a system's
planets, moons and belts, and the rows they leave in `system_interests`.

A reason is an entry of `REASONS`: an id (what the table, the wire format and the map's Menu call it), a test on `planets` (as `p`) and one
on `moons` (as `m`) that picks the bodies the reason counts, and optionally the column whose largest value over those bodies is kept as
the reason's score. A system has a row for each reason it has at least one body for. Adding a reason is one entry here (and one in
`static/galaxymarkers.js`'s `INTEREST_REASONS`); the table, the writers, the query and the node need no change.
"""

from planetgen import tuning
from planetgen.physics import habitability, life_classes

INHABITED = "inhabited"
CIVILIZATION = "civilization"
SENTIENT = "sentient"
LIFE = "life"
HABITABLE = "habitable"
UNINHABITABLE = "uninhabitable"
GAS_GIANT = "gas-giant"
ASTEROID_FIELD = "asteroid-field"

HABITABLE_MIN_SCORE = habitability.YELLOW_FLOOR
"""float: The lowest overall PHI-6 score of a world the map marks (Boss, 7:31 PM EDT 2026-10-24: "Yellow and above"). The habitable rule
itself stays Blue or Green (`habitability.is_habitable`, GEN.643); this is the map's reason, one tier wider. A scored terrestrial world
below it is an `UNINHABITABLE` one."""

TECH_REASONS = (CIVILIZATION,)
"""tuple: The reasons whose `best_score` is a technology level (0 to 7, `tuning.TECH_BANDS`) and not a habitability score: the map reads it as
the band of the system's most advanced species, not as a tier colour."""

BODY_TABLES = {"planets": ("planets", "p"), "moons": ("moons", "m"), "belts": ("asteroid_belts", "b")}
"""dict: The tables a reason's tests read, by the key a test is given under in `REASONS`: (table, alias)."""

# id -> {"planets"/"moons"/"belts": SQL test on `p`/`m`/`b` or None, "score": column of planets and moons or None, "params": values for
# the tests' `?`}. The reason's source, in one place each: a reason that reads a stored class instead (GEN.646) changes its entry here.
_CLASSES = ", ".join(f"'{letter}'" for letter in tuning.HABITABLE_PLANET_CLASSES)


def _inhabited_test(alias, host_column):
    """The test of a world that is inhabited the way the system page decides it (`query._with_life_fields`): its timeline reached a
    technological civilization (a habitable class only), or a colony stands on it (a terrestrial colony facility)."""
    return (f"(({alias}.life_stage = 'technological_civilization' AND {alias}.planet_class IN ({_CLASSES}))"
            f" OR EXISTS (SELECT 1 FROM facilities f WHERE f.{host_column} = {alias}.id AND f.kind = 'colony' AND f.placement = 'terrestrial'))")


def _tech_of(alias):
    """The score expression of the `civilization` reason (MAP.344): the technology level (`species.tech_level`) of the species whose homeworld
    is planet `alias`. A species lives on a planet, so a moon has none."""
    if alias != "p":
        return "NULL"
    return "(SELECT MAX(sp.tech_level) FROM species sp WHERE sp.homeworld_planet_id = p.id)"


def _life_test(alias, columns):
    """The test of a world with native life: at least one of the life scores `columns` (0 none to 9 abundant, NULL for a body with no score)
    is 1 or more."""
    return "(" + " OR ".join(f"{alias}.{column} >= 1" for column in columns) + ")"


_LIFE_COLUMNS = ("life_microbial", "life_flora", "life_fauna", "life_sentient")

REASONS = {
    INHABITED: {
        "planets": _inhabited_test("p", "planet_id"), "moons": _inhabited_test("m", "moon_id"), "score": "phi4", "params": (),
    },
    # Life (MAP.311): read from the stored life scores and the timeline's stage as they are, no expectation of how many worlds qualify. A
    # reason with no score keeps `best_score` NULL: the map's tier colour is the habitability score's, not a life score's. The civilization
    # reason's score is the technology level of its most advanced species (MAP.344, in `TECH_REASONS`), NULL for a world with none yet.
    CIVILIZATION: {
        "planets": f"p.life_stage = 'technological_civilization' AND p.planet_class IN ({_CLASSES})",
        "moons": f"m.life_stage = 'technological_civilization' AND m.planet_class IN ({_CLASSES})", "score": _tech_of, "params": (),
    },
    SENTIENT: {
        "planets": _life_test("p", ("life_sentient",)), "moons": _life_test("m", ("life_sentient",)), "score": None, "params": (),
    },
    LIFE: {
        "planets": _life_test("p", _LIFE_COLUMNS), "moons": _life_test("m", _LIFE_COLUMNS), "score": None, "params": (),
    },
    HABITABLE: {
        "planets": "p.phi4 >= ?", "moons": "m.phi4 >= ?", "score": "phi4", "params": (HABITABLE_MIN_SCORE,),
    },
    UNINHABITABLE: {
        # A scored terrestrial planet or moon under the habitable floor (MAP.321: the brown symbol).
        "planets": "p.body_type = 't' AND p.phi4 < ?", "moons": "m.phi4 < ?", "score": "phi4", "params": (HABITABLE_MIN_SCORE,),
    },
    GAS_GIANT: {
        # A planet whose class is a gas giant's (`body_type` 'g'); a moon is never one.
        "planets": "p.body_type = 'g'", "moons": None, "score": None, "params": (),
    },
    ASTEROID_FIELD: {
        # Any belt (MAP.321 rule 4; valuable ones with their resource rows later, MAP.312).
        "belts": "1 = 1", "score": None, "params": (),
    },
}
"""dict: The interest reasons, in the order the map prefers them (the first a system has names its mark, MAP.321: the best world's
symbol, an inhabited one before a civilization, a sentient species and any life, before a habitable one before a barren one before a gas giant
before a belt)."""

CLASS_PREFIX = "class-"
"""str: A planet class's reason is this and the class letter ("class-K"): a system with a planet or moon of that class (MAP.220, the
chosen classes). The visitor chooses which of them the map shows; they come after the reasons above, so a system that has one of those
keeps that mark."""


def class_reason(letter):
    """The reason id of planet class `letter`."""
    return CLASS_PREFIX + letter


for _letter in sorted(tuning.PLANET_CLASSES):
    REASONS[class_reason(_letter)] = {
        "planets": "p.planet_class IN (?)", "moons": "m.planet_class IN (?)", "score": None, "params": (_letter,),
    }
del _letter

LIFE_CLASS_PREFIX = "lifeclass-"
"""str: A life class's reason is this and the class id ("lifeclass-1a", "lifeclass-6"): a system with a planet or moon on which that life class
lives (`life_class` and `life_class_variant` on a body that has a `life_stage`, the class the `life` facet of the search reads, MAP.351, GEN.1001, or a row of its `planet_life` / `moon_life` for another class that lives there, GEN.1209). Like the planet classes, the visitor
chooses which of them the map shows, and they come after every reason above."""


def life_class_reason(class_id):
    """The reason id of life class `class_id` (such as `"1a"` or `"6"`)."""
    return LIFE_CLASS_PREFIX + class_id


for _definition in life_classes.CLASSES:
    REASONS[life_class_reason(_definition.id)] = {
        "planets": "((p.life_class = ? AND COALESCE(p.life_class_variant, '') = ? AND p.life_stage IS NOT NULL) "
                   "OR EXISTS (SELECT 1 FROM planet_life lc WHERE lc.planet_id = p.id AND lc.life_class = ? AND lc.life_stage IS NOT NULL))",
        "moons": "((m.life_class = ? AND COALESCE(m.life_class_variant, '') = ? AND m.life_stage IS NOT NULL) "
                 "OR EXISTS (SELECT 1 FROM moon_life lc WHERE lc.moon_id = m.id AND lc.life_class = ? AND lc.life_stage IS NOT NULL))", "score": None,
        "params": (_definition.number, _definition.variant or "", _definition.id),
    }
del _definition

CLASS_REASONS = tuple(reason for reason in REASONS if reason.startswith(CLASS_PREFIX))
"""tuple: The planet class reasons, one for each planet class, in letter order."""

LIFE_CLASS_REASONS = tuple(reason for reason in REASONS if reason.startswith(LIFE_CLASS_PREFIX))
"""tuple: The life class reasons, one for each class of `physics/life_classes.py`, in order."""

PER_CLASS_PREFIXES = (CLASS_PREFIX, LIFE_CLASS_PREFIX)
"""tuple: The prefixes of the reasons a node lists a few systems for each, and not by the bodies behind them (the planet and life classes)."""

BASE_REASONS = tuple(reason for reason in REASONS if not reason.startswith(PER_CLASS_PREFIXES))
"""tuple: The reasons that are neither a planet class nor a life class: what every node lists whatever the visitor chose."""


def reason_ids():
    """The reason ids, in order."""
    return tuple(REASONS)


def _sources(reason):
    """`[(key, test)]` of the tables a reason reads: its tests that are not None."""
    spec = REASONS[reason]
    return [(key, spec[key]) for key in BODY_TABLES if spec.get(key)]


def _params_of(reason):
    """The values for the `?` of a reason's tests, in the order `_count_sql` and `_score_sql` put the tests."""
    spec = REASONS[reason]
    return [value for _key, test in _sources(reason) for value in spec["params"] if test.count("?")]


def _count_sql(reason):
    """The SQL that counts a system's bodies for `reason` (a value for each `?` in `REASONS[reason]['params']`, tables in `BODY_TABLES` order)."""
    parts = []
    for key, test in _sources(reason):
        table, alias = BODY_TABLES[key]
        parts.append(f"(SELECT COUNT(*) FROM {table} {alias} WHERE {alias}.star_system_id = ss.id AND {test})")
    return " + ".join(parts), _params_of(reason)


def _score_sql(reason):
    """The SQL that gives a system's largest score among the bodies of `reason` (NULL for none or a reason with no score)."""
    spec = REASONS[reason]
    if not spec["score"]:
        return "NULL", []
    parts = []
    for key, test in _sources(reason):
        table, alias = BODY_TABLES[key]
        score = spec["score"](alias) if callable(spec["score"]) else f"{alias}.{spec['score']}"
        parts.append(f"(SELECT MAX({score}) FROM {table} {alias} WHERE {alias}.star_system_id = ss.id AND {test})")
    return "GREATEST(COALESCE(" + ", -1), COALESCE(".join(parts) + ", -1))", _params_of(reason)


def refresh_technology(conn):
    """
    Puts the `best_score` of the `TECH_REASONS` rows right (MAP.344): the technology level of the most advanced species on a world the reason
    counts. The species are drawn by the population pass after the planets are stored, and a species' level can change, so the pass and the
    migration that added the score call this; the rows themselves (which systems, how many worlds) do not depend on a species.
    """
    for reason in TECH_REASONS:
        score_sql, score_params = _score_sql(reason)
        conn.execute(
            "UPDATE system_interests i JOIN star_systems ss ON ss.id = i.star_system_id"
            f" SET i.best_score = NULLIF({score_sql}, -1) WHERE i.reason = ?", (*score_params, reason))


def record(conn, where, params):
    """
    Puts the `system_interests` rows of the systems matching `where` (SQL on `star_systems` aliased `ss`, with `params`) right: the rows
    they had go, and each system with a body for a reason gets a row for it with how many (at most 65535), the largest score and where
    the system is. A system with no galaxy position gets none.
    """
    conn.execute(
        "DELETE i FROM system_interests i JOIN star_systems ss ON ss.id = i.star_system_id"
        f" WHERE {where}", tuple(params))
    for reason in REASONS:
        count_sql, count_params = _count_sql(reason)
        score_sql, score_params = _score_sql(reason)
        conn.execute(
            "INSERT INTO system_interests (star_system_id, reason, n, best_score, galaxy_x_pc, galaxy_y_pc, galaxy_z_pc)"
            " SELECT id, ?, LEAST(65535, n), IF(best < 0, NULL, best), x, y, z FROM"
            f" (SELECT ss.id AS id, ss.galaxy_x_pc AS x, ss.galaxy_y_pc AS y, ss.galaxy_z_pc AS z, {count_sql} AS n, {score_sql} AS best"
            f" FROM star_systems ss WHERE ({where}) AND ss.galaxy_x_pc IS NOT NULL AND ss.galaxy_y_pc IS NOT NULL"
            " AND ss.galaxy_z_pc IS NOT NULL) c WHERE n > 0",
            (reason, *count_params, *score_params, *params))
