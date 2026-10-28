"""
The reasons a star system is worth a mark on the Galaxy Map (MAP.305, MAP.303): a registry of small SQL tests over a system's
planets and moons, and the rows they leave in `system_interests`.

A reason is an entry of `REASONS`: an id (what the table, the wire format and the map's Menu call it), a test on `planets` (as `p`) and one
on `moons` (as `m`) that picks the bodies the reason counts, and optionally the column whose largest value over those bodies is kept as
the reason's score. A system has a row for each reason it has at least one body for. Adding a reason is one entry here (and one in
`static/galaxymarkers.js`'s `INTEREST_REASONS`); the table, the writers, the query and the node need no change.
"""

from planetgen.physics import habitability

HABITABLE = "habitable"
GAS_GIANT = "gas-giant"

HABITABLE_MIN_SCORE = habitability.YELLOW_FLOOR
"""float: The lowest overall PHI-6 score of a world the map marks (Boss, 7:31 PM EDT 2026-10-24: "Yellow and above"). The habitable rule
itself stays Blue or Green (`habitability.is_habitable`, GEN.643); this is the map's reason, one tier wider."""

# id -> {"planets": SQL test on `p`, "moons": SQL test on `m` or None, "score": column or None, "params": values for the tests}
# The reason's source, in one place each: a reason that reads a stored class instead (GEN.646) changes its entry here.
REASONS = {
    HABITABLE: {
        "planets": "p.phi4 >= ?", "moons": "m.phi4 >= ?", "score": "phi4", "params": (HABITABLE_MIN_SCORE,),
    },
    GAS_GIANT: {
        # A planet whose class is a gas giant's (`body_type` 'g'); a moon is never one.
        "planets": "p.body_type = 'g'", "moons": None, "score": None, "params": (),
    },
}
"""dict: The interest reasons, in the order the map prefers them (the first a system has names its mark)."""


def reason_ids():
    """The reason ids, in order."""
    return tuple(REASONS)


def _count_sql(reason):
    """The SQL that counts a system's bodies for `reason` (a value for each `?` in `REASONS[reason]['params']`, planets then moons)."""
    spec = REASONS[reason]
    parts = [f"(SELECT COUNT(*) FROM planets p WHERE p.star_system_id = ss.id AND {spec['planets']})"]
    params = list(spec["params"]) if spec["planets"].count("?") else []
    if spec["moons"]:
        parts.append(f"(SELECT COUNT(*) FROM moons m WHERE m.star_system_id = ss.id AND {spec['moons']})")
        params += list(spec["params"]) if spec["moons"].count("?") else []
    return " + ".join(parts), params


def _score_sql(reason):
    """The SQL that gives a system's largest score among the bodies of `reason` (NULL for none or a reason with no score)."""
    spec = REASONS[reason]
    if not spec["score"]:
        return "NULL", []
    parts = [f"(SELECT MAX(p.{spec['score']}) FROM planets p WHERE p.star_system_id = ss.id AND {spec['planets']})"]
    params = list(spec["params"]) if spec["planets"].count("?") else []
    if spec["moons"]:
        parts.append(f"(SELECT MAX(m.{spec['score']}) FROM moons m WHERE m.star_system_id = ss.id AND {spec['moons']})")
        params += list(spec["params"]) if spec["moons"].count("?") else []
    return "GREATEST(COALESCE(" + ", -1), COALESCE(".join(parts) + ", -1))", params


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
