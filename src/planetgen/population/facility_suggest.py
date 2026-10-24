# planetgen/population/facility_suggest.py

"""
Suggested ratings for a facility (POP.10).

Nothing in the data holds a true crime figure or housing stock, so an admin sets ratings by hand. This is the
first rule set that fills the ones nobody set: each dimension gets a score from 0 to 1 out of what is stored
(the host's habitability, the six tech indices of the affiliated species, the kind of facility, the host's
richness, how far the system is from its polity's capital), a stable draw of up to `tuning.FACILITY_RATING_NOISE`
either way is added (seeded from the facility, so a rerun gives the same answer), and the score becomes Green at
`FACILITY_RATING_GREEN` or above, Yellow at `FACILITY_RATING_YELLOW` or above, Red below. See
`docs/design/habitability/population-and-politics.md`, "Setting ratings".

| Rating | Score |
|---|---|
| Health | 0.5 medical + 0.5 habitability |
| Housing | 0.5 habitability + 0.3 materials + 0.2 kind size |
| Resources | 0.4 host richness + 0.3 energy + 0.3 nearness |
| Maintenance | 0.35 materials + 0.25 energy + 0.25 nearness + 0.15 habitability |
| Crime | 0.5 information + 0.5 nearness, less a penalty when unaffiliated |

(An index counts as its value over 7; nearness is 1 less the system's distance from its polity's capital over
the polity's reach.)
"""

import hashlib

from planetgen import tuning
from planetgen.population import facility_ratings as ratings, tech

_WEIGHTS = {
    "health": {"medical": 0.5, "habitability": 0.5},
    "housing": {"habitability": 0.5, "materials": 0.3, "size": 0.2},
    "resources": {"richness": 0.4, "energy": 0.3, "nearness": 0.3},
    "maintenance": {"materials": 0.35, "energy": 0.25, "nearness": 0.25, "habitability": 0.15},
    "crime": {"information": 0.5, "nearness": 0.5},
}


def _noise(facility_id, dimension):
    """A stable draw between -`FACILITY_RATING_NOISE` and +`FACILITY_RATING_NOISE` for one rating of one facility."""
    digest = hashlib.sha256(f"facility-rating:{facility_id}:{dimension}".encode()).digest()
    unit = int.from_bytes(digest[:8], "big") / 2 ** 64
    return (unit * 2.0 - 1.0) * tuning.FACILITY_RATING_NOISE


def color_for(score):
    """The colour word of a 0 to 1 score."""
    if score >= tuning.FACILITY_RATING_GREEN:
        return "green"
    return "yellow" if score >= tuning.FACILITY_RATING_YELLOW else "red"


def scores(inputs, facility_id):
    """
    The five scores of a facility, 0 to 1, from `gather`'s inputs.

    Returns:
        dict: `{dimension: score}`, the noise included and the result kept inside 0 and 1.
    """
    found = {}
    for dimension, weights in _WEIGHTS.items():
        total = sum(weight * inputs[name] for name, weight in weights.items())
        if dimension == "crime" and not inputs["affiliated"]:
            total -= tuning.FACILITY_RATING_UNAFFILIATED_PENALTY
        found[dimension] = min(1.0, max(0.0, total + _noise(facility_id, dimension)))
    return found


def suggest(inputs, facility_id):
    """`{dimension: colour word}` for a facility's `gather`ed inputs."""
    return {dimension: color_for(score) for dimension, score in scores(inputs, facility_id).items()}


def gather(conn, facility_id):
    """
    What the rules read about one facility, each input 0 to 1 (`tuning.FACILITY_RATING_UNKNOWN` when it is not
    stored): `habitability`, `richness`, `size`, `nearness`, the six tech indices of the affiliated species by name
    (`energy`, `materials`, `information`, `medical`, `propulsion`, `defense`) and `affiliated`. `None` for an
    unknown facility.
    """
    unknown = tuning.FACILITY_RATING_UNKNOWN
    row = conn.execute(
        "SELECT kind, placement, host_type, planet_id, moon_id, star_system_id, affiliation_polity_id FROM facilities "
        "WHERE id = ?", (facility_id,)).fetchone()
    if row is None:
        return None
    habitability = tuning.FACILITY_RATING_ORBITAL_HABITABILITY
    if row["placement"] == "terrestrial":
        table, key = ("planets", row["planet_id"]) if row["planet_id"] else ("moons", row["moon_id"])
        stored = conn.execute(f"SELECT phi4 FROM {table} WHERE id = ?", (key,)).fetchone()
        habitability = unknown if stored is None or stored["phi4"] is None else float(stored["phi4"])
    inputs = {
        "habitability": habitability,
        "richness": tuning.FACILITY_HOST_RICHNESS.get(row["host_type"], unknown),
        "size": tuning.FACILITY_KIND_SIZE.get(row["kind"], unknown),
        "nearness": 0.0,
        "affiliated": row["affiliation_polity_id"] is not None,
        **{domain: unknown for domain in tech.DOMAINS},
    }
    polity_id = row["affiliation_polity_id"]
    if polity_id is not None:
        polity = conn.execute(
            "SELECT p.reach_ly, " + ", ".join(f"s.{domain}_index" for domain in tech.DOMAINS) + " FROM polities p "
            "JOIN species s ON s.id = p.species_id WHERE p.id = ?", (polity_id,)).fetchone()
        if polity is not None:
            for domain in tech.DOMAINS:
                index = polity[f"{domain}_index"]
                inputs[domain] = unknown if index is None else index / tech.MAX_INDEX
            owner = conn.execute("SELECT polity_id, distance_ly FROM system_owners WHERE star_system_id = ?",
                                 (row["star_system_id"],)).fetchone() if row["star_system_id"] else None
            if owner is not None and owner["polity_id"] == polity_id and polity["reach_ly"]:
                inputs["nearness"] = min(1.0, max(0.0, 1.0 - owner["distance_ly"] / polity["reach_ly"]))
    return inputs


def suggest_for_facility(conn, facility_id):
    """
    Fills the unset ratings of a facility with suggestions.

    Only a rating that is NULL changes. The facility's `ratings_source` becomes `generated` when it had no
    ratings before or they were all generated; with any set by hand it stays `admin`.

    Returns:
        list: The dimensions filled (empty when none was unset), or `None` for an unknown facility.
    """
    inputs = gather(conn, facility_id)
    if inputs is None:
        return None
    row = conn.execute("SELECT " + ", ".join(f"rating_{name}" for name in ratings.RATINGS) + ", ratings_source "
                       "FROM facilities WHERE id = ?", (facility_id,)).fetchone()
    unset = [name for name in ratings.RATINGS if row[f"rating_{name}"] is None]
    if not unset:
        return []
    chosen = suggest(inputs, facility_id)
    columns = {f"rating_{name}": ratings.value_of(chosen[name]) for name in unset}
    clauses = [f"`{column}` = ?" for column in columns]
    if len(unset) == len(ratings.RATINGS) or row["ratings_source"] == "generated":
        clauses.append("ratings_source = 'generated'")
    clauses.append("ratings_set_at = CURRENT_TIMESTAMP")
    conn.execute(f"UPDATE facilities SET {', '.join(clauses)} WHERE id = ?", (*columns.values(), facility_id))
    return unset


BATCH_LIMIT = 5000
"""int: The most facilities one batch of suggestions fills."""


def suggest_for_unassessed(conn):
    """
    Fills the unset ratings of every facility that has any, up to `BATCH_LIMIT` facilities.

    Returns:
        tuple: `(facilities changed, ratings filled)`.
    """
    ids = [row["id"] for row in conn.execute(
        "SELECT id FROM facilities WHERE " + " OR ".join(f"rating_{name} IS NULL" for name in ratings.RATINGS)
        + " ORDER BY id LIMIT ?", (BATCH_LIMIT,)).fetchall()]
    facilities = filled = 0
    for facility_id in ids:
        done = suggest_for_facility(conn, facility_id)
        if done:
            facilities += 1
            filled += len(done)
    return facilities, filled
