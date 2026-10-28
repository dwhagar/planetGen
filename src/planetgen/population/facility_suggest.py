# planetgen/population/facility_suggest.py

"""
Suggested ratings for a facility (POP.10).

Nothing in the data holds a true crime figure or housing stock, so an admin sets ratings by hand. This is the
first rule set that fills the ones nobody set: each dimension gets a score from 0 to 1 out of what is stored
(the host's habitability, Suitability, radiation domain, food rung and Flora and Fauna scores, the six tech indices of
the affiliated species, the kind of facility, the host's richness, how far the system is from its polity's capital), a
stable draw of up to `tuning.FACILITY_RATING_NOISE`
either way is added (seeded from the facility, so a rerun gives the same answer), and the score becomes Green at
`FACILITY_RATING_GREEN` or above, Yellow at `FACILITY_RATING_YELLOW` or above, Red below. See
`docs/design/habitability/population-and-politics.md`, "Setting ratings".

| Rating | Score |
|---|---|
| Health | 0.4 medical + 0.2 habitability + 0.2 radiation + 0.2 food |
| Housing | 0.25 habitability + 0.25 Suitability + 0.3 materials + 0.2 kind size |
| Resources | 0.4 host richness + 0.3 energy + 0.3 nearness; on a living world 0.3 host richness + 0.2 biota + 0.25 energy + 0.25 nearness |
| Maintenance | 0.35 materials + 0.25 energy + 0.25 nearness + 0.15 habitability |
| Crime | 0.5 information + 0.5 nearness, less a penalty when unaffiliated |

(An index counts as its value over 7; nearness is 1 less the system's distance from its polity's capital over
the polity's reach. Habitability is the host's PHI-4, Suitability and radiation are its `phi4_suitability` and
`phi4_radiation`, food is the Suitability food rung worked out from its stored life scores, and biota is its Flora
and Fauna scores over 18; a host that is not a world, a world without a score and a lifeless world have no biota, and
the lifeless one keeps the first Resources weights. An orbital host reads `FACILITY_RATING_ORBITAL_HABITABILITY` for
every host condition.)
"""

import hashlib

from planetgen import tuning
from planetgen.db import storage_checksum, store
from planetgen.physics import habitability, life_scores
from planetgen.population import facility_ratings as ratings, tech

_WEIGHTS = {
    "health": {"medical": 0.4, "habitability": 0.2, "radiation": 0.2, "food": 0.2},
    "housing": {"habitability": 0.25, "suitability": 0.25, "materials": 0.3, "size": 0.2},
    "resources": {"richness": 0.4, "energy": 0.3, "nearness": 0.3},
    "maintenance": {"materials": 0.35, "energy": 0.25, "nearness": 0.25, "habitability": 0.15},
    "crime": {"information": 0.5, "nearness": 0.5},
}

_LIVING_RESOURCES = {"richness": 0.3, "biota": 0.2, "energy": 0.25, "nearness": 0.25}
"""dict: The Resources weights on a living world (`gather` gives it a `biota`): its Flora and Fauna are a renewable resource."""

MAX_LIFE_SCORE = 9
"""int: The top of a Flora or Fauna score (`life_scores`), so biota is their sum over twice this."""


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
        if dimension == "resources" and inputs.get("biota") is not None:
            weights = _LIVING_RESOURCES
        total = sum(weight * inputs[name] for name, weight in weights.items())
        if dimension == "crime" and not inputs["affiliated"]:
            total -= tuning.FACILITY_RATING_UNAFFILIATED_PENALTY
        found[dimension] = min(1.0, max(0.0, total + _noise(facility_id, dimension)))
    return found


def suggest(inputs, facility_id):
    """`{dimension: colour word}` for a facility's `gather`ed inputs."""
    return {dimension: color_for(score) for dimension, score in scores(inputs, facility_id).items()}


def _score(value, unknown):
    """A stored 0 to 1 score, or `unknown` when the body has none."""
    return unknown if value is None else float(value)


def _food(conn, table, stored, unknown):
    """The Suitability food rung (`habitability.food_score`) of a stored world, from the stage its Microbial, Flora and
    Fauna scores give (microbes under ice or in a deep aquifer are out of reach), or `unknown` for a body with no score."""
    world = store.stored_world(conn, table, stored, with_life=False)
    if world is None:
        return unknown
    world.life_stage = life_scores.food_stage(stored["life_microbial"], stored["life_flora"], stored["life_fauna"])
    return habitability.food_score(world)


def gather(conn, facility_id):
    """
    What the rules read about one facility, each input 0 to 1 (`tuning.FACILITY_RATING_UNKNOWN` when it is not
    stored): `habitability`, `suitability`, `radiation`, `food`, `biota` (`None` unless the host is a living world),
    `richness`, `size`, `nearness`, the six tech indices of the affiliated species by name (`energy`, `materials`,
    `information`, `medical`, `propulsion`, `defense`) and `affiliated`. `None` for an unknown facility.
    """
    unknown = tuning.FACILITY_RATING_UNKNOWN
    row = conn.execute(
        "SELECT kind, placement, host_type, planet_id, moon_id, star_system_id, affiliation_polity_id FROM facilities "
        "WHERE id = ?", (facility_id,)).fetchone()
    if row is None:
        return None
    host = dict.fromkeys(("habitability", "suitability", "radiation", "food"), tuning.FACILITY_RATING_ORBITAL_HABITABILITY)
    biota = None
    if row["placement"] == "terrestrial":
        table, key = ("planets", row["planet_id"]) if row["planet_id"] else ("moons", row["moon_id"])
        stored = conn.execute(f"SELECT * FROM {table} WHERE id = ?", (key,)).fetchone()
        host = dict.fromkeys(host, unknown)
        if stored is not None and stored["phi4"] is not None:
            host.update(habitability=float(stored["phi4"]), suitability=_score(stored["phi4_suitability"], unknown),
                        radiation=_score(stored["phi4_radiation"], unknown), food=_food(conn, table, stored, unknown))
            flora, fauna = stored["life_flora"] or 0, stored["life_fauna"] or 0
            if flora or fauna:
                biota = (flora + fauna) / (2 * MAX_LIFE_SCORE)
    inputs = {
        **host,
        "biota": biota,
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
    storage_checksum.mark_rows(conn, "facilities", [facility_id])
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
