# planetgen/population/model.py

"""
Population and politics (POP.1 to POP.4, schema v44)
====================================================

Names the species of every world with a technological civilization
(GEN.80; worlds with simpler life get none), dates each one and places it in an era, founds a polity for every
spacefaring species, and splits the generated systems around them into
territories. See `docs/design/population-and-politics.md`.

Everything here works from what is already stored -- the planets'
evolutionary paragraphs (`evolution.get_evolutionary_timeline`) and the
systems' positions -- so it runs as its own pass (`planetgen
population`, and after every `galaxy`/`sector` run) on a galaxy that was
generated before it existed. Nothing in system generation changes.

The pure functions (`parse_timeline`, `species_traits`,
`civilization_age`, `era_for_age`, `reach_ly`, `resolve_claims`, ...)
take no database; `run_pass` and the functions after it do.
"""

import colorsys
import contextlib
import hashlib
import math
import re
from dataclasses import dataclass

import pymysql

from planetgen import tuning
from planetgen.util import draw
from planetgen.util import log
from planetgen.generation.evolution import life_stage_from_paragraphs
from planetgen.names.wordlists import STAR_NAMES, STAR_PREFIXES, STAR_SUFFIXES
from planetgen.names.wordsalad import generate_phoneme_salad_name
from planetgen.physics.units import ly_to_pc, pc_to_ly

_YEARS_PER_UNIT = {"Billion": 1e9, "Million": 1e6}

_SYSTEM_AGE_RE = re.compile(r"current estimated age of the system is ([0-9.]+) (Billion|Million) Years")
_MILESTONE_AGE_RE = re.compile(r"would have been [A-Za-z ]+ at ([0-9.]+) (Billion|Million) Years")

SCAN_BATCH = 5000
"""int: Planet ids per scan query."""

_NAME_ATTEMPTS = 50
"""int: Fresh draws before a species name collision gives up and takes a
numbered suffix instead."""

PASS_LOCK_WAIT_SECONDS = 10
"""int: How long one wait for another population pass's lock lasts
before `run_pass` says it is waiting (and waits again)."""


# ---------------------------------------------------------------------------
# Pure model
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class Timeline:
    """What a planet's evolutionary paragraphs say: the milestone it
    reached and the system age and milestone age, in years."""
    life_stage: str
    system_age_years: float
    milestone_age_years: float

    @property
    def window_years(self):
        """Years since the milestone: how old a civilization could be."""
        return max(0.0, self.system_age_years - self.milestone_age_years)


def _years(match):
    return float(match.group(1)) * _YEARS_PER_UNIT[match.group(2)]


def parse_timeline(paragraphs):
    """
    Reads a planet's stored evolutionary paragraphs.

    Args:
        paragraphs (list[str]): `planet_evolutionary_paragraphs` rows, in
            order.

    Returns:
        Timeline or None: `None` unless the milestone is one of
            `tuning.LIFE_WORLD_STAGES` and both ages can be
            read back.
    """
    stage = life_stage_from_paragraphs(paragraphs)
    if stage not in tuning.LIFE_WORLD_STAGES:
        return None
    text = " ".join(paragraphs)
    system_age = _SYSTEM_AGE_RE.search(text)
    milestone_age = _MILESTONE_AGE_RE.search(text)
    if system_age is None or milestone_age is None:
        return None
    return Timeline(stage, _years(system_age), _years(milestone_age))


def species_traits(gravity_g, surface_temperature_k, rng):
    """
    A dominant species' build, climate and size, from its homeworld.

    Args:
        gravity_g (float or None): Surface gravity in g.
        surface_temperature_k (float or None): Surface temperature in K.
        rng (draw.Stream): The draw for size.

    Returns:
        dict: `build`, `climate` and `size`.
    """
    gravity = gravity_g if gravity_g is not None else 1.0
    temperature = surface_temperature_k if surface_temperature_k is not None else 288.0
    build = "gracile" if gravity < 0.7 else "robust" if gravity > 1.4 else "medium"
    climate = "cold" if temperature < 250 else "hot" if temperature > 320 else "temperate"
    # Heavier worlds favor smaller bodies, lighter worlds larger ones.
    weights = {"gracile": (1, 2, 3), "medium": (1, 2, 1), "robust": (3, 2, 1)}[build]
    size = rng.choices(("small", "medium", "large"), weights=weights, k=1)[0]
    return {"build": build, "climate": climate, "size": size}


def has_civilization(timeline, forced, rng, prevalence_percent=0.0):
    """
    Whether a life world has a technological civilization now: its
    timeline must have reached the milestone, and then either intelligent
    life was forced on for its system (`forced`, `system_configs.
    intelligent_life`) or a `CIVILIZATION_CHANCE` draw succeeds, that
    chance moved by the generating run's intelligent-life prevalence
    (`prevalence_percent`, `system_configs.prevalence_intelligent_life`,
    GEN.52).
    """
    if timeline.life_stage != "technological_civilization":
        return False
    if forced:
        return True
    chance = min(1.0, max(0.0, tuning.CIVILIZATION_CHANCE * (1.0 + (prevalence_percent or 0.0) / 100.0)))
    return rng.random() < chance


def civilization_age(window_years, rng):
    """
    A technological civilization's age in years: log-uniform between
    `CIVILIZATION_MIN_AGE_YEARS` and `window_years`, so ages spread evenly
    across orders of magnitude and never exceed what the star allows. A
    window shorter than the minimum (the stored ages are rounded to 10
    million years, so this is rounding) gives the minimum.
    """
    low = tuning.CIVILIZATION_MIN_AGE_YEARS
    if window_years <= low:
        return low
    return math.exp(rng.uniform(math.log(low), math.log(window_years)))


def era_for_age(age_years):
    """`(era, spacefaring)` for a civilization `age_years` old
    (`tuning.CIVILIZATION_ERAS`)."""
    era, spacefaring = None, False
    for name, start, is_spacefaring in tuning.CIVILIZATION_ERAS:
        if age_years >= start:
            era, spacefaring = name, is_spacefaring
    return era, spacefaring


def interstellar_start_years():
    """The age at which a civilization first counts as spacefaring."""
    return min(start for _name, start, spacefaring in tuning.CIVILIZATION_ERAS if spacefaring)


def reach_ly(age_years, cap_ly=None):
    """
    A polity's territory reach: `TERRITORY_BASE_REACH_LY` when it goes
    interstellar, growing with the square root of age, capped at `cap_ly`
    (`TERRITORY_REACH_CAP_LY` by default). 0 before it is spacefaring.
    """
    cap = tuning.TERRITORY_REACH_CAP_LY if cap_ly is None else cap_ly
    start = interstellar_start_years()
    if age_years < start:
        return 0.0
    return min(cap, tuning.TERRITORY_BASE_REACH_LY * math.sqrt(age_years / start))


def government_for(species_id):
    """The polity's form of government, stable for a species."""
    forms = tuning.GOVERNMENT_FORMS
    return forms[draw.Stream(species_id * 7919).randrange(len(forms))]


def polity_name(species_name, government):
    """e.g. `"Voranthis Concord"`."""
    return f"{species_name} {government}"


def polity_color(species_id):
    """A map color for a polity (`#rrggbb`): golden-ratio hues so
    neighbouring ids differ, saturated enough to read on both themes."""
    hue = (species_id * 0.618033988749895) % 1.0
    red, green, blue = colorsys.hls_to_rgb(hue, 0.55, 0.7)
    return "#{:02x}{:02x}{:02x}".format(round(red * 255), round(green * 255), round(blue * 255))


def new_species_name():
    """A fresh species name candidate, from the star-name generator."""
    return generate_phoneme_salad_name(STAR_NAMES, STAR_PREFIXES, STAR_SUFFIXES, allow_split=False, max_length=10)


def add_claims(best, polity_id, capital, reach, systems):
    """
    Adds one polity's claims to `best` (`system_id -> (strength,
    polity_id, distance_ly)`), keeping the stronger claim on each system.
    A claim's strength is `reach / distance`; a capital's own system (at
    distance 0) is infinitely strong. Ties keep the claim already there,
    so callers go in ascending polity id.

    Args:
        best (dict): The running claims, updated in place.
        polity_id (int): The claimant.
        capital (tuple): `(x, y, z)` pc.
        reach (float): Light years.
        systems (dict): `system_id -> (x, y, z)` pc.
    """
    if reach <= 0:
        return
    reach_pc = ly_to_pc(reach)
    for system_id, position in systems.items():
        distance_pc = math.dist(capital, position)
        if distance_pc > reach_pc:
            continue
        strength = math.inf if distance_pc == 0 else reach_pc / distance_pc
        held = best.get(system_id)
        if held is None or strength > held[0]:
            best[system_id] = (strength, polity_id, pc_to_ly(distance_pc))


def resolve_claims(polities, systems):
    """
    Splits systems among polities: each system inside some polity's
    reach goes to the strongest claim (`add_claims`). Ties go to the lower
    polity id.

    Args:
        polities (list[tuple]): `(polity_id, capital (x, y, z) pc, reach_ly)`.
        systems (dict): `system_id -> (x, y, z)` pc, every candidate.

    Returns:
        dict: `system_id -> (polity_id, distance_ly)`.
    """
    best = {}
    for polity_id, capital, reach in sorted(polities):
        add_claims(best, polity_id, capital, reach, systems)
    return {system_id: (polity_id, distance) for system_id, (_s, polity_id, distance) in best.items()}


# ---------------------------------------------------------------------------
# Database pass
# ---------------------------------------------------------------------------

def run_pass(conn, rescan=False, territories_only=False):
    """
    The whole population pass, committed as it goes: scan new life
    worlds, refresh eras and polities, recompute territories.

    Args:
        conn (Connection): A write connection.
        rescan (bool): Forget every species (and with them every polity
            and territory) and scan every planet again.
        territories_only (bool): Skip straight to recomputing territories.

    Returns:
        dict: Counts -- `new_species`, `species`, `spacefaring`,
            `polities`, `owned_systems`.
    """
    with _one_pass_at_a_time(conn):
        return _run_pass(conn, rescan, territories_only)


def _pass_lock_name(conn):
    """The named lock for this database's population pass (named locks
    are server-wide, so the database is part of the name)."""
    database = conn.execute("SELECT DATABASE() AS db").fetchone()["db"] or ""
    return "planetgen-population:" + hashlib.sha256(database.encode("utf-8")).hexdigest()[:32]


@contextlib.contextmanager
def _one_pass_at_a_time(conn):
    """
    Holds this database's population lock (`GET_LOCK`) for a whole pass,
    so two passes at once -- a `--population` run finishing while
    another, or `planetgen population`, is still going (TEST.37) --
    take turns: the second waits, then sees the first's species and
    watermark, rather than both naming the same worlds and one failing
    on a duplicate homeworld or species name.
    """
    name = _pass_lock_name(conn)
    waited = False
    while True:
        got = conn.execute("SELECT GET_LOCK(?, ?) AS ok", (name, PASS_LOCK_WAIT_SECONDS)).fetchone()["ok"]
        if got == 1:
            break
        if got is None:
            raise pymysql.err.OperationalError(1205, f"could not take the population lock {name!r}")
        if not waited:
            waited = True
            log.normal("Population: waiting for another population pass to finish.")
    # A fresh transaction, so the pass reads what the last holder committed.
    conn.commit()
    try:
        yield
    finally:
        try:
            conn.execute("SELECT RELEASE_LOCK(?)", (name,))
        except pymysql.err.MySQLError:  # a lost session has already dropped it
            pass


def _run_pass(conn, rescan, territories_only):
    new_species = 0
    if not territories_only:
        if rescan:
            with conn:
                conn.execute("DELETE FROM species")
                conn.execute("DELETE FROM population_state")
        new_species = scan_life_worlds(conn)
        drop_species_without_civilization(conn)
        refresh_civilizations(conn)
    owned = refresh_territories(conn)
    counts = conn.execute(
        "SELECT COUNT(*) AS species, COALESCE(SUM(spacefaring), 0) AS spacefaring FROM species"
    ).fetchone()
    polities = conn.execute("SELECT COUNT(*) AS n FROM polities").fetchone()["n"]
    return {"new_species": new_species, "species": int(counts["species"]),
            "spacefaring": int(counts["spacefaring"]), "polities": int(polities), "owned_systems": owned}


def _watermark(conn):
    row = conn.execute("SELECT scanned_planet_id FROM population_state WHERE id = 1").fetchone()
    return row["scanned_planet_id"] if row else 0


def _set_watermark(conn, planet_id):
    conn.execute(
        "INSERT INTO population_state (id, scanned_planet_id) VALUES (1, ?) "
        "ON DUPLICATE KEY UPDATE scanned_planet_id = VALUES(scanned_planet_id)", (planet_id,),
    )


def _unique_species_name(conn, taken):
    """A species name no stored species (or `taken`) has."""
    for _attempt in range(_NAME_ATTEMPTS):
        name = new_species_name()
        if name in taken:
            continue
        if conn.execute("SELECT 1 FROM species WHERE name = ?", (name,)).fetchone() is None:
            return name
    base = new_species_name()
    for number in range(2, 10_000):
        name = f"{base} {number}"
        if name not in taken and conn.execute("SELECT 1 FROM species WHERE name = ?", (name,)).fetchone() is None:
            return name
    raise RuntimeError("could not find a free species name")


def scan_life_worlds(conn):
    """
    Adds a `species` row for every world with a technological civilization
    (GEN.80) among the planets not yet scanned
    (`population_state.scanned_planet_id`), batch by batch, and moves the
    watermark past them.

    Returns:
        int: Species added.
    """
    added = 0
    top = conn.execute("SELECT MAX(id) AS id FROM planets").fetchone()["id"] or 0
    low = _watermark(conn)
    while low < top:
        high = min(top, low + SCAN_BATCH)
        rows = conn.execute(
            "SELECT p.id, p.star_system_id, p.life_chemical, p.gravity_g, p.surface_temperature_k, "
            "cfg.intelligent_life, cfg.prevalence_intelligent_life, e.paragraph FROM planets p "
            "JOIN planet_evolutionary_paragraphs e ON e.planet_id = p.id "
            "JOIN star_systems ss ON ss.id = p.star_system_id "
            "JOIN system_configs cfg ON cfg.id = ss.system_config_id "
            "LEFT JOIN species s ON s.homeworld_planet_id = p.id "
            "WHERE p.id > ? AND p.id <= ? AND s.id IS NULL "
            "AND e.paragraph LIKE ? ORDER BY p.id, e.position",
            (low, high, "%would have been Technological Civilization at%"),
        ).fetchall()
        planets = {}
        for row in rows:
            planets.setdefault(row["id"], (row, []))[1].append(row["paragraph"])
        taken = set()
        inserts = []
        for planet_id, (row, paragraphs) in planets.items():
            timeline = parse_timeline(paragraphs)
            if timeline is None:
                continue
            rng = draw.Stream(planet_id)
            traits = species_traits(row["gravity_g"], row["surface_temperature_k"], rng)
            if not has_civilization(timeline, row["intelligent_life"], rng, row["prevalence_intelligent_life"]):
                # Only a technological civilization gets a species (GEN.80).
                continue
            age = civilization_age(timeline.window_years, rng)
            era, spacefaring = era_for_age(age)
            name = _unique_species_name(conn, taken)
            taken.add(name)
            inserts.append((name, planet_id, row["star_system_id"], row["life_chemical"], timeline.life_stage,
                            traits["build"], traits["climate"], traits["size"], age, era, int(spacefaring)))
        with conn:
            if inserts:
                conn.executemany(
                    "INSERT INTO species (name, homeworld_planet_id, star_system_id, life_chemical, life_stage, "
                    "build, climate, size, civilization_age_years, era, spacefaring) "
                    "VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)", inserts,
                )
            _set_watermark(conn, high)
        added += len(inserts)
        low = high
    if added:
        log.debug(f"Population: {added} new species")
    return added


def drop_species_without_civilization(conn):
    """
    Deletes the species an earlier pass named for worlds with no
    technological civilization (GEN.80: only those get a species). They
    never had a polity; anything else of theirs goes with them (`ON DELETE
    CASCADE`).

    Returns:
        int: Species deleted.
    """
    with conn:
        dropped = conn.execute("DELETE FROM species WHERE civilization_age_years IS NULL").rowcount
    if dropped:
        log.debug(f"Population: {dropped} species without a technological civilization removed")
    return dropped


def refresh_civilizations(conn):
    """
    Recomputes every civilization's era from its stored age (picking up
    retuned `CIVILIZATION_ERAS`), founds a polity for each newly
    spacefaring species, updates every polity's reach, and dissolves
    polities whose species is no longer spacefaring.
    """
    with conn:
        rows = conn.execute(
            "SELECT id, civilization_age_years, era, spacefaring FROM species "
            "WHERE civilization_age_years IS NOT NULL"
        ).fetchall()
        updates = []
        for row in rows:
            era, spacefaring = era_for_age(row["civilization_age_years"])
            if era != row["era"] or int(spacefaring) != row["spacefaring"]:
                updates.append((era, int(spacefaring), row["id"]))
        if updates:
            conn.executemany("UPDATE species SET era = ?, spacefaring = ? WHERE id = ?", updates)
        conn.execute(
            "DELETE p FROM polities p JOIN species s ON s.id = p.species_id WHERE s.spacefaring = 0"
        )
        founding = conn.execute(
            "SELECT s.id, s.name, s.star_system_id FROM species s "
            "LEFT JOIN polities p ON p.species_id = s.id WHERE s.spacefaring = 1 AND p.id IS NULL"
        ).fetchall()
        if founding:
            conn.executemany(
                "INSERT INTO polities (name, species_id, capital_system_id, government, color, reach_ly) "
                "VALUES (?, ?, ?, ?, ?, 0)",
                [(polity_name(row["name"], government_for(row["id"])), row["id"], row["star_system_id"],
                  government_for(row["id"]), polity_color(row["id"])) for row in founding],
            )
        reaches = conn.execute(
            "SELECT p.id, p.reach_ly, s.civilization_age_years FROM polities p JOIN species s ON s.id = p.species_id"
        ).fetchall()
        changed = [(reach_ly(row["civilization_age_years"]), row["id"]) for row in reaches
                   if not math.isclose(reach_ly(row["civilization_age_years"]), row["reach_ly"])]
        if changed:
            conn.executemany("UPDATE polities SET reach_ly = ? WHERE id = ?", changed)


def _system_positions(conn, sector_ids):
    """`system_id -> (x, y, z)` galaxy-frame parsecs for every placed
    system in `sector_ids`."""
    positions = {}
    ids = sorted(sector_ids)
    for start in range(0, len(ids), 500):
        chunk = ids[start:start + 500]
        marks = ", ".join("?" * len(chunk))
        rows = conn.execute(
            "SELECT ss.id, ss.position_x_mpc + sec.center_x_pc * 1000 AS x, "
            "ss.position_y_mpc + sec.center_y_pc * 1000 AS y, "
            "ss.position_z_mpc + sec.center_z_pc * 1000 AS z "
            "FROM star_systems ss JOIN sectors sec ON sec.id = ss.sector_id "
            f"WHERE ss.sector_id IN ({marks}) AND ss.position_x_mpc IS NOT NULL "
            "AND sec.center_x_pc IS NOT NULL",
            tuple(chunk),
        ).fetchall()
        for row in rows:
            positions[row["id"]] = (row["x"] / 1000.0, row["y"] / 1000.0, row["z"] / 1000.0)
    return positions


def capital_positions(conn):
    """`[(polity_id, (x, y, z) pc or None, reach_ly)]` for every polity."""
    rows = conn.execute(
        "SELECT p.id, p.reach_ly, ss.position_x_mpc + sec.center_x_pc * 1000 AS x, "
        "ss.position_y_mpc + sec.center_y_pc * 1000 AS y, ss.position_z_mpc + sec.center_z_pc * 1000 AS z "
        "FROM polities p JOIN star_systems ss ON ss.id = p.capital_system_id "
        "LEFT JOIN sectors sec ON sec.id = ss.sector_id ORDER BY p.id"
    ).fetchall()
    return [(row["id"], None if row["x"] is None else (row["x"] / 1000.0, row["y"] / 1000.0, row["z"] / 1000.0),
             row["reach_ly"]) for row in rows]


def refresh_territories(conn):
    """
    Recomputes `system_owners` from scratch. A polity whose capital was
    never placed in the galaxy (a standalone system, or a sector made
    with no galaxy) holds only its capital.

    Returns:
        int: Systems owned.
    """
    from planetgen.db.store import sectors_reached_by

    # Each polity only looks at the sectors its own reach touches, so the
    # cost is the sum of their territories, not polities x all systems.
    best = {}
    for polity_id, capital, reach in capital_positions(conn):
        if capital is None or reach <= 0:
            continue
        sectors = sectors_reached_by(conn, capital, ly_to_pc(reach))
        add_claims(best, polity_id, capital, reach, _system_positions(conn, sectors))
    owners = {system_id: (polity_id, distance) for system_id, (_s, polity_id, distance) in best.items()}
    # Every capital holds its own system, placed or not.
    for row in conn.execute(
        "SELECT p.id, p.capital_system_id FROM polities p ORDER BY p.id DESC"
    ).fetchall():
        held = owners.get(row["capital_system_id"])
        if held is None or held[1] > 0:
            owners[row["capital_system_id"]] = (row["id"], 0.0)
    with conn:
        conn.execute("DELETE FROM system_owners")
        rows = [(system_id, polity_id, distance) for system_id, (polity_id, distance) in owners.items()]
        for start in range(0, len(rows), 5000):
            conn.executemany(
                "INSERT INTO system_owners (star_system_id, polity_id, distance_ly) VALUES (?, ?, ?)",
                rows[start:start + 5000],
            )
    return len(owners)


# ---------------------------------------------------------------------------
# Reads (API, pages)
# ---------------------------------------------------------------------------

def population_status(conn):
    """
    What population data exists, for pages to hide themselves when there
    is none (Boss, 2026-10-01): `generated` (a pass has run), `species`
    (any species), `polities` (any polity) and `territories` (any owned
    system). Three indexed `EXISTS` probes. All false on a database
    without the v44 tables.
    """
    try:
        row = conn.execute(
            "SELECT EXISTS(SELECT 1 FROM population_state) AS `generated`, "
            "EXISTS(SELECT 1 FROM species) AS species, EXISTS(SELECT 1 FROM polities) AS polities, "
            "EXISTS(SELECT 1 FROM system_owners) AS territories"
        ).fetchone()
    except Exception:  # noqa: BLE001 -- a pre-v44 database: nothing generated
        return {"generated": False, "species": False, "polities": False, "territories": False}
    return {key: bool(row[key]) for key in ("generated", "species", "polities", "territories")}


_SPECIES_COLUMNS = (
    "s.id, s.name, s.homeworld_planet_id, p.name AS homeworld_name, s.star_system_id, ss.name AS system_name, "
    "s.life_chemical, s.life_stage, s.build, s.climate, s.size, s.civilization_age_years, s.era, s.spacefaring, "
    "pol.id AS polity_id, pol.name AS polity_name"
)
_SPECIES_FROM = (
    "FROM species s JOIN planets p ON p.id = s.homeworld_planet_id "
    "JOIN star_systems ss ON ss.id = s.star_system_id LEFT JOIN polities pol ON pol.species_id = s.id"
)


def _species_dict(row):
    item = dict(row)
    item["spacefaring"] = bool(item["spacefaring"])
    return item


SPECIES_SORTS = {
    "name": "s.name", "homeworld": "p.name", "system": "ss.name", "era": "s.era", "spacefaring": "s.spacefaring",
    "polity": "pol.name",
}
"""dict: The sort keys `list_species` accepts (the Species table's column keys) -> the SQL they order by."""


def _species_where(spacefaring=None, eras=()):
    """`(where_sql, params)` for the Species table's filters (an empty `eras` filters nothing)."""
    clauses, params = [], []
    if spacefaring is not None:
        clauses.append("s.spacefaring = ?")
        params.append(int(spacefaring))
    if eras:
        clauses.append(f"s.era IN ({', '.join('?' for _ in eras)})")
        params.extend(eras)
    return ("WHERE " + " AND ".join(clauses) + " ") if clauses else "", params


def list_species(conn, spacefaring=None, limit=50, offset=0, sort="name", descending=False, eras=()):
    """A page of species, by `sort` (a key of `SPECIES_SORTS`; ties fall back
    to name, then id), optionally only (non-)spacefaring ones or those in
    one of the `eras`. A species with no polity sorts last by `"polity"`."""
    if sort not in SPECIES_SORTS:
        raise ValueError(f"unknown species sort {sort!r}")
    where, params = _species_where(spacefaring, eras)
    column = SPECIES_SORTS[sort]
    nulls_last = f"{column} IS NULL, " if sort in ("polity", "era") else ""
    rows = conn.execute(
        f"SELECT {_SPECIES_COLUMNS} {_SPECIES_FROM} {where}"
        f"ORDER BY {nulls_last}{column} {'DESC' if descending else 'ASC'}, s.name, s.id LIMIT ? OFFSET ?",
        params + [limit, offset],
    ).fetchall()
    return [_species_dict(row) for row in rows]


def count_species(conn, spacefaring=None, eras=()):
    where, params = _species_where(spacefaring, eras)
    return conn.execute(f"SELECT COUNT(*) AS n FROM species s {where}", params).fetchone()["n"]


def species_facets(conn, spacefaring=None, eras=()):
    """
    The option counts for the Species table's menus: `spacefaring` (`"yes"`/
    `"no"`) and `era`. Each menu ignores its own filter.

    Returns:
        dict: `{"spacefaring"|"era": [{"value", "count"}]}`, options with no species left out.
    """
    def counts(select, **filters):
        where, params = _species_where(**{"spacefaring": spacefaring, "eras": eras, **filters})
        return {row["value"]: row["n"] for row in conn.execute(
            f"SELECT {select} AS value, COUNT(*) AS n FROM species s {where}GROUP BY value", params).fetchall()}

    flags = counts("(CASE WHEN s.spacefaring THEN 'yes' ELSE 'no' END)", spacefaring=None)
    eras_found = counts("s.era", eras=())
    return {
        "spacefaring": [{"value": v, "count": flags[v]} for v in ("yes", "no") if flags.get(v)],
        "era": [{"value": v, "count": eras_found[v]} for v in sorted(x for x in eras_found if x is not None)],
    }


def species_detail(conn, species_id):
    """One species, or `None`."""
    row = conn.execute(f"SELECT {_SPECIES_COLUMNS} {_SPECIES_FROM} WHERE s.id = ?", (species_id,)).fetchone()
    return None if row is None else _species_dict(row)


def species_on_planet(conn, planet_id):
    """The species whose homeworld is this planet, or `None`."""
    row = conn.execute(
        f"SELECT {_SPECIES_COLUMNS} {_SPECIES_FROM} WHERE s.homeworld_planet_id = ?", (planet_id,)
    ).fetchone()
    return None if row is None else _species_dict(row)


_POLITY_COLUMNS = (
    "pol.id, pol.name, pol.government, pol.color, pol.reach_ly, pol.species_id, s.name AS species_name, "
    "s.era, s.civilization_age_years, pol.capital_system_id, ss.name AS capital_name, "
    "(SELECT COUNT(*) FROM system_owners o WHERE o.polity_id = pol.id) AS system_count"
)
_POLITY_FROM = (
    "FROM polities pol JOIN species s ON s.id = pol.species_id JOIN star_systems ss ON ss.id = pol.capital_system_id"
)


POLITY_SORTS = {
    "name": "pol.name", "species": "s.name", "government": "pol.government", "capital": "ss.name",
    "systems": "system_count", "reach": "pol.reach_ly",
}
"""dict: The sort keys `list_polities` accepts (the Polities table's column keys) -> the SQL they order by."""


def _polity_where(governments=(), eras=()):
    """`(where_sql, params)` for the Polities table's filters (the species' `eras`)."""
    clauses, params = [], []
    for column, values in (("pol.government", governments), ("s.era", eras)):
        if values:
            clauses.append(f"{column} IN ({', '.join('?' for _ in values)})")
            params.extend(values)
    return ("WHERE " + " AND ".join(clauses) + " ") if clauses else "", params


def list_polities(conn, limit=50, offset=0, sort="name", descending=False, governments=(), eras=()):
    """A page of polities, by `sort` (a key of `POLITY_SORTS`; ties fall back
    to name, then id), optionally only those with one of the `governments`
    or whose species is in one of the `eras`."""
    if sort not in POLITY_SORTS:
        raise ValueError(f"unknown polity sort {sort!r}")
    where, params = _polity_where(governments, eras)
    rows = conn.execute(
        f"SELECT {_POLITY_COLUMNS} {_POLITY_FROM} {where}"
        f"ORDER BY {POLITY_SORTS[sort]} {'DESC' if descending else 'ASC'}, pol.name, pol.id LIMIT ? OFFSET ?",
        params + [limit, offset],
    ).fetchall()
    return [dict(row) for row in rows]


def count_polities(conn, governments=(), eras=()):
    where, params = _polity_where(governments, eras)
    return conn.execute(
        f"SELECT COUNT(*) AS n FROM polities pol JOIN species s ON s.id = pol.species_id {where}", params
    ).fetchone()["n"]


def polity_facets(conn, governments=(), eras=()):
    """
    The option counts for the Polities table's menus, `government` and `era`
    (the species'). Each menu ignores its own filter.

    Returns:
        dict: `{"government"|"era": [{"value", "count"}]}`, options with no polities left out.
    """
    def counts(select, **filters):
        where, params = _polity_where(**{"governments": governments, "eras": eras, **filters})
        return {row["value"]: row["n"] for row in conn.execute(
            f"SELECT {select} AS value, COUNT(*) AS n FROM polities pol "
            f"JOIN species s ON s.id = pol.species_id {where}GROUP BY value", params).fetchall()}

    found_governments = counts("pol.government", governments=())
    found_eras = counts("s.era", eras=())
    return {
        "government": [{"value": v, "count": found_governments[v]}
                       for v in sorted(x for x in found_governments if x is not None)],
        "era": [{"value": v, "count": found_eras[v]} for v in sorted(x for x in found_eras if x is not None)],
    }


POLITY_SYSTEM_SORTS = {"name": "ss.name", "distance": "o.distance_ly"}
"""dict: The sort keys `polity_detail` accepts for its systems."""


def polity_detail(conn, polity_id, limit=50, offset=0, sort="distance", descending=False):
    """One polity with a page of the systems it owns (nearest the capital
    first unless sorted by a key of `POLITY_SYSTEM_SORTS`), or `None`."""
    if sort not in POLITY_SYSTEM_SORTS:
        raise ValueError(f"unknown polity system sort {sort!r}")
    row = conn.execute(f"SELECT {_POLITY_COLUMNS} {_POLITY_FROM} WHERE pol.id = ?", (polity_id,)).fetchone()
    if row is None:
        return None
    systems = conn.execute(
        "SELECT o.star_system_id AS id, ss.name, o.distance_ly FROM system_owners o "
        "JOIN star_systems ss ON ss.id = o.star_system_id WHERE o.polity_id = ? "
        f"ORDER BY {POLITY_SYSTEM_SORTS[sort]} {'DESC' if descending else 'ASC'}, o.distance_ly, "
        "o.star_system_id LIMIT ? OFFSET ?",
        (polity_id, limit, offset),
    ).fetchall()
    return {**dict(row), "systems": [dict(system) for system in systems], "limit": limit, "offset": offset}


def system_owner(conn, star_system_id):
    """`{polity_id, polity_name, color, distance_ly}` for the polity that
    owns a system, or `None`."""
    row = conn.execute(
        "SELECT o.polity_id, pol.name AS polity_name, pol.color, o.distance_ly FROM system_owners o "
        "JOIN polities pol ON pol.id = o.polity_id WHERE o.star_system_id = ?", (star_system_id,),
    ).fetchone()
    return None if row is None else dict(row)


def territory_points(conn, limit=20000):
    """
    Owned systems with galaxy-frame positions (parsecs) and polity colors,
    for a 3D territory overlay -- at most `limit` rows, capitals first.
    """
    rows = conn.execute(
        "SELECT o.star_system_id AS id, o.polity_id, pol.color, "
        "(ss.position_x_mpc + sec.center_x_pc * 1000) / 1000 AS x, "
        "(ss.position_y_mpc + sec.center_y_pc * 1000) / 1000 AS y, "
        "(ss.position_z_mpc + sec.center_z_pc * 1000) / 1000 AS z "
        "FROM system_owners o JOIN polities pol ON pol.id = o.polity_id "
        "JOIN star_systems ss ON ss.id = o.star_system_id JOIN sectors sec ON sec.id = ss.sector_id "
        "WHERE ss.position_x_mpc IS NOT NULL AND sec.center_x_pc IS NOT NULL "
        "ORDER BY o.distance_ly, o.star_system_id LIMIT ?", (limit,),
    ).fetchall()
    return [dict(row) for row in rows]
