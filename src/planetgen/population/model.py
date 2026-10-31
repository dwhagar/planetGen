# planetgen/population/model.py

"""
Population and politics (POP.1 to POP.4, schema v44)
====================================================

Names the species of every world with a technological civilization
(GEN.80; worlds with simpler life get none), dates each one and places it in an era, founds a polity for every
spacefaring species, and splits the generated systems around them into
territories. See `docs/design/habitability/population-and-politics.md`.

Everything here works from what is already stored -- the planets'
evolutionary paragraphs (`evolution.get_evolutionary_state`) and the
systems' positions -- so it runs as its own pass (`planetgen
population` and `planetgen territories`, and as stages of a `galaxy` run) on a galaxy that was
generated before it existed. Nothing in system generation changes.

The pure functions (`timeline_of`, `species_traits`,
`civilization_age`, `era_for_age`, `reach_ly`, `resolve_claims`, ...)
take no database; `run_population`, `run_territories` and the functions after them do.
"""

import colorsys
import contextlib
import hashlib
import math
from dataclasses import dataclass

import pymysql

from planetgen import tuning
from planetgen.population import tech
from planetgen.util import draw
from planetgen.util import generator_epoch
from planetgen.util import log, request_odds
from planetgen.names.wordlists import STAR_NAMES, STAR_PREFIXES, STAR_SUFFIXES
from planetgen.names.wordsalad import generate_phoneme_salad_name
from planetgen.physics import life_classes
from planetgen.physics.units import ly_to_pc, pc_to_ly
from planetgen.util import detmath
from planetgen.physics import unit_registry as units

_mpc_to_pc = units.converter("mpc", "pc")

SCAN_BATCH = 5000
"""int: Planet ids per scan query."""

_NAME_ATTEMPTS = 50
"""int: Fresh draws before a species name collision gives up and takes a
numbered suffix instead."""

PASS_LOCK_WAIT_SECONDS = 10
"""int: How long one wait for another population pass's lock lasts
before a pass says it is waiting (and waits again)."""


# ---------------------------------------------------------------------------
# Pure model
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class Timeline:
    """What a planet's stored timeline says: the milestone it
    reached and the system age and milestone age, in years."""
    life_stage: str
    system_age_years: float
    milestone_age_years: float

    @property
    def window_years(self):
        """Years since the milestone: how old a civilization could be."""
        return max(0.0, self.system_age_years - self.milestone_age_years)


def timeline_of(life_stage, system_age_gyr, milestone_age_gyr):
    """
    A planet's stored timeline (`planets.life_stage`, `life_system_age_gyr`, `life_milestone_age_gyr`).

    Returns:
        Timeline or None: `None` unless the milestone is one of `tuning.LIFE_WORLD_STAGES` and both ages are
            stored.
    """
    if life_stage not in tuning.LIFE_WORLD_STAGES or system_age_gyr is None or milestone_age_gyr is None:
        return None
    return Timeline(life_stage, system_age_gyr * 1e9, milestone_age_gyr * 1e9)


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


def has_civilization(timeline, requested, rng, prevalence_percent=0.0):
    """
    Whether a life world has a technological civilization now: its
    timeline must have reached the milestone, and then a `CIVILIZATION_CHANCE`
    draw succeeds, that chance moved by the generating run's intelligent-life
    prevalence (`prevalence_percent`, `system_configs.prevalence_intelligent_life`,
    GEN.52). A request for intelligent life on its system (`requested`,
    `system_configs.intelligent_life`) multiplies the odds of that draw by
    `tuning.REQUEST_ODDS` and forces nothing (GEN.1130).
    """
    if timeline.life_stage != "technological_civilization":
        return False
    chance = min(1.0, max(0.0, tuning.CIVILIZATION_CHANCE * (1.0 + (prevalence_percent or 0.0) / 100.0)))
    if requested:
        chance = request_odds.raised_chance(chance)
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


def object_key(uid):
    """
    The text a population stream is keyed on for a world: its unique ID in hex (GEN.306). Row ids come from the
    id blocks of however many workers filled the galaxy, so a stream keyed on one would draw differently at 1 and
    at 4 workers; the unique ID is given once at birth and is the same however the galaxy was filled.
    """
    if uid is None:
        raise ValueError("a world with no unique ID cannot be given a stream of its own")
    return bytes(uid).hex()


def government_for(homeworld_uid):
    """The polity's form of government, stable for a species: drawn from its homeworld's unique ID."""
    forms = tuning.GOVERNMENT_FORMS
    return forms[draw.Stream(f"government:{object_key(homeworld_uid)}").randrange(len(forms))]


def polity_name(species_name, government):
    """e.g. `"Voranthis Concord"`."""
    return f"{species_name} {government}"


def polity_color(homeworld_uid):
    """A map color for a polity (`#rrggbb`): a hue drawn from its homeworld's unique ID, saturated enough to
    read on both themes."""
    hue = draw.Stream(f"polity-color:{object_key(homeworld_uid)}").random()
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
        distance_pc = detmath.dist(capital, position)
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

def run_population(conn, rescan=False, on_stage=None):
    """
    The population pass, committed as it goes (PERF.125): scan new life worlds, drop the species without a
    civilization, date the civilizations and found or dissolve the polities. It does not draw territories; that is
    `run_territories`, its own command.

    Args:
        conn (Connection): A write connection.
        rescan (bool): Forget every species (and with them every polity
            and territory) and scan every planet again.
        on_stage (callable, optional): `on_stage(description)` as each of
            the pass's `POPULATION_STAGES` steps starts (UX.83's bar).

    Returns:
        dict: Counts -- `new_species`, `species`, `spacefaring`, `polities`.
    """
    with _one_pass_at_a_time(conn):
        return _run_population(conn, rescan, on_stage)


def run_territories(conn, on_stage=None):
    """
    The territories pass, committed as it goes (PERF.125): recompute which polity owns which system from the polities
    already stored.

    Args:
        conn (Connection): A write connection.
        on_stage (callable, optional): `on_stage(description)` as the one step starts.

    Returns:
        dict: Counts -- `polities` and `owned_systems`.
    """
    with _one_pass_at_a_time(conn):
        if on_stage is not None:
            on_stage("Population: recomputing territories")
        owned = refresh_territories(conn)
        polities = conn.execute("SELECT COUNT(*) AS n FROM polities").fetchone()["n"]
        return {"polities": int(polities), "owned_systems": owned}


def _pass_lock_name(conn):
    """The named lock for this database's population pass (named locks
    are server-wide, so the database is part of the name)."""
    database = conn.execute("SELECT DATABASE() AS db").fetchone()["db"] or ""
    return "planetgen-population:" + hashlib.sha256(database.encode("utf-8")).hexdigest()[:32]


@contextlib.contextmanager
def _one_pass_at_a_time(conn):
    """
    Holds this database's population lock (`GET_LOCK`) for a whole pass,
    so two passes at once -- a `galaxy` run's population stage finishing while
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


POPULATION_STAGES = 3
"""int: How many steps `run_population` announces through `on_stage`."""


def _run_population(conn, rescan, on_stage=None):
    def stage(description):
        if on_stage is not None:
            on_stage(description)

    stage("Population: scanning life worlds")
    if rescan:
        with conn:
            conn.execute("DELETE FROM species")
            conn.execute("DELETE FROM population_state")
    new_species = scan_life_worlds(conn)
    stage("Population: dropping species without a civilization")
    drop_species_without_civilization(conn)
    stage("Population: dating civilizations and founding polities")
    refresh_civilizations(conn)
    from planetgen.db import interests  # here: the db package imports this one

    with conn:
        interests.refresh_technology(conn)    # MAP.344: the map's civilization rows carry their most advanced species' level
    counts = conn.execute(
        "SELECT COUNT(*) AS species, COALESCE(SUM(spacefaring), 0) AS spacefaring FROM species"
    ).fetchone()
    polities = conn.execute("SELECT COUNT(*) AS n FROM polities").fetchone()["n"]
    return {"new_species": new_species, "species": int(counts["species"]),
            "spacefaring": int(counts["spacefaring"]), "polities": int(polities)}


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
    from planetgen.db import species_fit  # here: the db package imports this one

    added = 0
    top = conn.execute("SELECT MAX(id) AS id FROM planets").fetchone()["id"] or 0
    low = _watermark(conn)
    while low < top:
        high = min(top, low + SCAN_BATCH)
        rows = conn.execute(
            "SELECT p.id, p.uid, p.star_system_id, p.life_class, p.life_class_variant, p.gravity_g, p.surface_temperature_k, "
            "p.life_stage, p.life_system_age_gyr, p.life_milestone_age_gyr, "
            "cfg.intelligent_life, cfg.prevalence_intelligent_life FROM planets p "
            "JOIN star_systems ss ON ss.id = p.star_system_id "
            "JOIN system_configs cfg ON cfg.id = ss.system_config_id "
            "LEFT JOIN species s ON s.homeworld_planet_id = p.id "
            "WHERE p.id > ? AND p.id <= ? AND s.id IS NULL "
            "AND p.life_stage = 'technological_civilization' ORDER BY p.uid",
            (low, high),
        ).fetchall()
        taken = set()
        inserts = []
        for row in rows:
            planet_id = row["id"]
            timeline = timeline_of(row["life_stage"], row["life_system_age_gyr"], row["life_milestone_age_gyr"])
            if timeline is None:
                continue
            key = object_key(row["uid"])
            rng = draw.Stream(f"life:{key}")
            traits = species_traits(row["gravity_g"], row["surface_temperature_k"], rng)
            if not has_civilization(timeline, row["intelligent_life"], rng, row["prevalence_intelligent_life"]):
                # Only a technological civilization gets a species (GEN.80).
                continue
            try:
                species_fit.check_homeworld(conn, planet_id)  # POP.14: only where its own limits let it live
            except species_fit.PlacementError as error:
                log.debug(f"Population: no species on planet {planet_id}: {error}")
                continue
            age = civilization_age(timeline.window_years, rng)
            era, spacefaring = era_for_age(age)
            with draw.bound(f"species-name:{key}"):
                name = _unique_species_name(conn, taken)
            taken.add(name)
            inserts.append((name, planet_id, row["star_system_id"], life_classes.class_id(row["life_class"], row["life_class_variant"]), timeline.life_stage,
                            traits["build"], traits["climate"], traits["size"], age, era, int(spacefaring),
                            generator_epoch.generator_epoch()))
        with conn:
            if inserts:
                conn.executemany(
                    "INSERT INTO species (name, homeworld_planet_id, star_system_id, life_class, life_stage, "
                    "build, climate, size, civilization_age_years, era, spacefaring, generator_epoch) "
                    "VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)", inserts,
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


TECH_INDEX_COLUMNS = tuple(f"{domain}_index" for domain in tech.DOMAINS)
"""tuple: The six `species` index columns, in `tuning.TECH_DOMAINS` order."""


def refresh_tech_levels(conn):
    """
    Draws the six tech indices and the tech level of every species that has none (POP.9): a new species, or
    one that existed before the columns. A species whose `tech_source` is `'admin'` (edited by hand) is never
    redrawn. The draw uses a stream of its own (`tech:<homeworld unique ID>`), apart from the one
    `scan_life_worlds` draws a species' traits and age from, so no stored age changes.

    Returns:
        int: Species drawn.
    """
    missing = " OR ".join(f"{column} IS NULL" for column in TECH_INDEX_COLUMNS)
    rows = conn.execute(
        "SELECT s.id, pl.uid, s.civilization_age_years FROM species s "
        "JOIN planets pl ON pl.id = s.homeworld_planet_id "
        "WHERE s.civilization_age_years IS NOT NULL AND (s.tech_source IS NULL OR s.tech_source <> 'admin') "
        f"AND ({missing})"
    ).fetchall()
    updates = []
    for row in rows:
        indices = tech.draw_indices(row["civilization_age_years"], draw.Stream(f"tech:{object_key(row['uid'])}"))
        updates.append((*(indices[domain] for domain in tech.DOMAINS), tech.tech_level(indices), row["id"]))
    if updates:
        sets = ", ".join(f"{column} = ?" for column in TECH_INDEX_COLUMNS)
        with conn:
            conn.executemany(
                f"UPDATE species SET {sets}, tech_level = ?, tech_source = 'generated' WHERE id = ?", updates)
        log.debug(f"Population: {len(updates)} species given tech levels")
    return len(updates)


def refresh_civilizations(conn):
    """
    Recomputes every civilization's era from its stored age (picking up
    retuned `CIVILIZATION_ERAS`), founds a polity for each newly
    spacefaring species, updates every polity's reach, and dissolves
    polities whose species is no longer spacefaring. Species with no tech indices get theirs
    (`refresh_tech_levels`).
    """
    refresh_tech_levels(conn)
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
            "SELECT s.id, s.name, s.star_system_id, pl.uid FROM species s "
            "JOIN planets pl ON pl.id = s.homeworld_planet_id "
            "LEFT JOIN polities p ON p.species_id = s.id WHERE s.spacefaring = 1 AND p.id IS NULL"
        ).fetchall()
        if founding:
            conn.executemany(
                "INSERT INTO polities (name, species_id, capital_system_id, government, color, reach_ly) "
                "VALUES (?, ?, ?, ?, ?, 0)",
                [(polity_name(row["name"], government_for(row["uid"])), row["id"], row["star_system_id"],
                  government_for(row["uid"]), polity_color(row["uid"])) for row in founding],
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
            positions[row["id"]] = (_mpc_to_pc(row["x"]), _mpc_to_pc(row["y"]), _mpc_to_pc(row["z"]))
    return positions


def capital_positions(conn):
    """`[(polity_id, (x, y, z) pc or None, reach_ly)]` for every polity."""
    rows = conn.execute(
        "SELECT p.id, p.reach_ly, ss.position_x_mpc + sec.center_x_pc * 1000 AS x, "
        "ss.position_y_mpc + sec.center_y_pc * 1000 AS y, ss.position_z_mpc + sec.center_z_pc * 1000 AS z "
        "FROM polities p JOIN star_systems ss ON ss.id = p.capital_system_id "
        "LEFT JOIN sectors sec ON sec.id = ss.sector_id ORDER BY p.id"
    ).fetchall()
    return [(row["id"], None if row["x"] is None else (_mpc_to_pc(row["x"]), _mpc_to_pc(row["y"]), _mpc_to_pc(row["z"])),
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
    "s.life_class, s.life_stage, s.build, s.climate, s.size, s.civilization_age_years, s.era, s.spacefaring, "
    "pol.id AS polity_id, pol.name AS polity_name, "
    "s.energy_index, s.materials_index, s.information_index, s.medical_index, s.propulsion_index, s.defense_index, "
    "s.tech_level, s.tech_source"
)
_SPECIES_FROM = (
    "FROM species s JOIN planets p ON p.id = s.homeworld_planet_id "
    "JOIN star_systems ss ON ss.id = s.star_system_id LEFT JOIN polities pol ON pol.species_id = s.id"
)


def _species_dict(row):
    item = dict(row)
    item["spacefaring"] = bool(item["spacefaring"])
    indices = {domain: item[f"{domain}_index"] for domain in tech.DOMAINS}
    if item["tech_level"] is None or None in indices.values():
        item["tech_band"] = item["limiting_domain"] = item["lopsided"] = None
    else:
        item["tech_level"] = float(item["tech_level"])
        item["tech_band"] = tech.band_name(item["tech_level"])
        item["limiting_domain"] = tech.limiting_domain(indices)
        item["lopsided"] = tech.is_lopsided(indices)
    return item


SPECIES_SORTS = {
    "name": "s.name", "homeworld": "p.name", "system": "ss.name", "era": "s.era", "spacefaring": "s.spacefaring",
    "polity": "pol.name", "tech": "s.tech_level",
}
"""dict: The sort keys `list_species` accepts (the Species table's column keys) -> the SQL they order by."""


TECH_RANGE_COLUMNS = ("tech_level",) + TECH_INDEX_COLUMNS
"""tuple: The columns `list_species` can search by range (`tech`)."""


def _species_where(spacefaring=None, eras=(), tech_ranges=None):
    """`(where_sql, params)` for the Species table's filters (an empty `eras` filters nothing). `tech_ranges`
    maps a `TECH_RANGE_COLUMNS` name to `(low, high)`, either end `None` for open."""
    clauses, params = [], []
    for column, (low, high) in (tech_ranges or {}).items():
        if column not in TECH_RANGE_COLUMNS:
            raise ValueError(f"unknown tech range {column!r}")
        if low is not None:
            clauses.append(f"s.{column} >= ?")
            params.append(low)
        if high is not None:
            clauses.append(f"s.{column} <= ?")
            params.append(high)
    if spacefaring is not None:
        clauses.append("s.spacefaring = ?")
        params.append(int(spacefaring))
    if eras:
        clauses.append(f"s.era IN ({', '.join('?' for _ in eras)})")
        params.extend(eras)
    return ("WHERE " + " AND ".join(clauses) + " ") if clauses else "", params


def _and(where, *clauses):
    """`where` (empty or `WHERE ... ` with its trailing space) with more conditions ANDed on."""
    clauses = [clause for clause in clauses if clause]
    if not clauses:
        return where
    return (where + "AND " if where else "WHERE ") + " AND ".join(clauses) + " "


SPECIES_WALK_MIN_ROWS = 50_000
"""int: Under this many species (the highest id) the `homeworld` and `system` sorts sort what they match, which costs next to nothing; from it up
they walk the planets or the systems by name (`_species_walk`, PERF.211)."""

SPECIES_CHUNK_MIN = 200
SPECIES_CHUNK_MAX = 20_000

_SPECIES_WALKS = {"homeworld": ("planets", "p", "s.homeworld_planet_id"), "system": ("star_systems", "ss", "s.star_system_id")}
"""dict: The sorts that order by a name in another table -> the table walked by name, its alias in `_SPECIES_FROM`, and the species column that
names a row of it."""


def _species_walk(conn, limit, sort, descending, where, params, after):
    """
    The rows of a page of species sorted by the name of its homeworld or its system, read in index order (PERF.211): the planets or
    systems are read by name a chunk at a time, only the species of those are sorted (a planet has one species at most), and the chunk
    doubles while the page is not full, so a page costs about what it returns whatever its depth. Ties go by the id of the planet or
    system, then the species' name and id, all in the direction of the sort. `None` when `after` is no species (the caller starts at
    the top the slow way).
    """
    from planetgen.db import query      # not at the top: the query module reads this package's other modules
    table, alias, link = _SPECIES_WALKS[sort]
    direction = "DESC" if descending else "ASC"
    start = None
    if after is not None:
        start = conn.execute(f"SELECT {alias}.name AS owner_name, {alias}.id AS owner_id, s.name AS name {_SPECIES_FROM} WHERE s.id = ?",
                             (after,)).fetchone()
        if start is None:
            return None
    rows, beyond, chunk = [], None, min(max(limit * 2, SPECIES_CHUNK_MIN), SPECIES_CHUNK_MAX)
    while len(rows) < limit:
        if beyond is not None:
            condition, owners_params = query._lexicographic(["name", "id"], list(beyond), descending)
            owners_where = f"WHERE {condition} "
        elif start is not None:
            condition, owners_params = query._lexicographic(["name", "id"], [start["owner_name"], start["owner_id"]], descending,
                                                            inclusive=True)
            owners_where = f"WHERE {condition} "
        else:
            owners_where, owners_params = "", []
        owners = conn.execute(f"SELECT id, name FROM {table} {owners_where}ORDER BY name {direction}, id {direction} LIMIT ?",
                              owners_params + [chunk]).fetchall()
        if not owners:
            break
        marks = ", ".join("?" for _ in owners)
        page_where = _and(where, f"{link} IN ({marks})")
        page_params = list(params) + [owner["id"] for owner in owners]
        if start is not None and beyond is None:
            condition, found = query._lexicographic([f"{alias}.name", f"{alias}.id", "s.name", "s.id"],
                                                    [start["owner_name"], start["owner_id"], start["name"], after], descending)
            page_where = _and(page_where, condition)
            page_params += found
        order = ", ".join(f"{key} {direction}" for key in (f"{alias}.name", f"{alias}.id", "s.name", "s.id"))
        rows += conn.execute(f"SELECT {_SPECIES_COLUMNS} {_SPECIES_FROM} {page_where}ORDER BY {order} LIMIT ?",
                             page_params + [limit - len(rows)]).fetchall()
        if len(owners) < chunk:
            break
        beyond = (owners[-1]["name"], owners[-1]["id"])
        chunk = min(chunk * 2, SPECIES_CHUNK_MAX)
    return rows


def list_species(conn, spacefaring=None, limit=50, offset=0, sort="name", descending=False, eras=(),
                 tech_ranges=None, after=None):
    """A page of species, by `sort` (a key of `SPECIES_SORTS`; ties fall back to name, in the same direction), optionally
    only (non-)spacefaring ones or those in one of the `eras`. A species with no polity, era or tech level sorts last by
    that sort.

    `after` is the id of the species the previous page ended on: the page starts after it in the same order and `offset`
    is not used, so a deep page costs what the first does (PERF.205). The `name`, `era`, `tech` and `spacefaring` sorts
    are read in index order (`idx_species_*_name`), the species with a value and then the ones with none, as two walks;
    the `polity` sort walks the few species that have a polity, and the `homeworld` and `system` sorts (a name in another table) walk the planets or
    the systems by name and sort the species of a few at a time (PERF.211; a page by `offset` past the first, or with fewer than
    `SPECIES_WALK_MIN_ROWS` species, sorts what it matches). An id that is not a species starts at the top.
    """
    if sort not in SPECIES_SORTS:
        raise ValueError(f"unknown species sort {sort!r}")
    from planetgen.db import query     # not at the top: the query module reads this package's other modules
    where, params = _species_where(spacefaring, eras, tech_ranges)
    column = SPECIES_SORTS[sort]
    if sort in _SPECIES_WALKS and (after is not None or not offset) and \
            conn.execute("SELECT COALESCE(MAX(id), 0) AS n FROM species").fetchone()["n"] >= SPECIES_WALK_MIN_ROWS:
        walked = _species_walk(conn, limit, sort, descending, where, params, after)
        if walked is not None:
            return [_species_dict(row) for row in walked]
    nullable = sort in ("polity", "era", "tech")
    keys = ["s.name"] if sort == "name" else [column, "s.name"]
    walks = [(f"{column} IS NOT NULL", keys), (f"{column} IS NULL", ["s.name"])] if nullable else [("", keys)]
    direction = "DESC" if descending else "ASC"
    first_walk, start = 0, None
    if after is not None:
        found = conn.execute(f"SELECT {column} AS value, s.name AS name {_SPECIES_FROM} WHERE s.id = ?", (after,)).fetchone()
        if found is not None:
            first_walk = 1 if nullable and found["value"] is None else 0
            start = found
            offset = 0
    rows = []
    for number, (clause, walk_keys) in enumerate(walks):
        if number < first_walk:
            continue
        walk_where = _and(where, clause)
        walk_params = list(params)
        if start is not None and number == first_walk:
            values = [start["name"]] if len(walk_keys) == 1 else [start["value"], start["name"]]
            condition, found_params = query._keyset_condition(walk_keys, values, descending, "s.id", after)
            walk_where = _and(walk_where, condition)
            walk_params += found_params
        if offset and number < len(walks) - 1:
            size = conn.execute(f"SELECT COUNT(*) AS n {_SPECIES_FROM} {walk_where}", walk_params).fetchone()["n"]
            if offset >= size:
                offset -= size
                continue
        order = ", ".join(f"{key} {direction}" for key in walk_keys) + f", s.id {direction}"
        room = limit - len(rows)
        rows += conn.execute(
            f"SELECT {_SPECIES_COLUMNS} {_SPECIES_FROM} {walk_where}ORDER BY {order} LIMIT ? OFFSET ?",
            walk_params + [room, offset]).fetchall()
        offset = 0
        if len(rows) >= limit:
            break
    return [_species_dict(row) for row in rows]


def count_species(conn, spacefaring=None, eras=(), tech_ranges=None):
    where, params = _species_where(spacefaring, eras, tech_ranges)
    return conn.execute(f"SELECT COUNT(*) AS n FROM species s {where}", params).fetchone()["n"]


def species_facets(conn, spacefaring=None, eras=(), tech_ranges=None):
    """
    The option counts for the Species table's menus: `spacefaring` (`"yes"`/
    `"no"`) and `era`. Each menu ignores its own filter.

    Returns:
        dict: `{"spacefaring"|"era": [{"value", "count"}]}`, options with no species left out.
    """
    def counts(select, **filters):
        where, params = _species_where(**{"spacefaring": spacefaring, "eras": eras, "tech_ranges": tech_ranges,
                                          **filters})
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
