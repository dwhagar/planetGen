# planetgen/galaxy/wide_binary_reroll.py

"""
Re-rolling the wide binaries that break the stability rule (GEN.229).

GEN.228 made new wide (S-type) binaries stable: each star's planets stay inside its stable orbit limit, and the two stars' outermost bodies keep clear of each other at the pair's
closest approach (`validation.check_star_system`). A galaxy made before that holds wide binaries whose orbital paths overlap and would end in planetary collisions. Boss decided
(8:37 AM EDT 2026-10-24): re-roll only the affected systems from their seeds, not repair them in place, and print to the console which systems were redone.

`unstable` reads every wide binary and returns those that fail the rule. `reroll` makes one of them again: the same kind of system (a wide binary, drawn for the same sector's stellar
population and distance from the galactic center), from a seed of its own, `unit_seed(galaxy seed, "system-reroll", <the system's ID>)`, so the same galaxy gives the same new system on
every run and a re-rolled system, which now follows the rule, is never touched again. The system keeps its row, name, ID, place and links (`store.replace_star_system_content`); its
bodies are made again under the slots of the old ones (a body in the same place keeps its ID). `run` does the lot, with a line for each system, and is the stage the bring-to-epoch run
and `planetgen bring-to-epoch --wide-binaries` share.

A system is left alone, and named, when it hosts facilities (they would be deleted with its bodies), is built around a bright star the galaxy scatter placed, or sits in a sector an
admin edited (the bring-to-epoch run's own rule, GEN.764).
"""

from planetgen.db import store
from planetgen.galaxy import mass_chain, object_uid, seed as galaxy_seed
from planetgen.generation import sector_context, validation
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics import constants

SEED_KIND = "system-reroll"
"""str: The kind of the unit seed a re-rolled system is drawn from (`galaxy.seed.unit_seed`); it must never change, or a galaxy would re-roll to a different system."""

PROGRESS_EVERY = 1000
"""int: The scan says where it is after this many wide binaries."""

MAX_ATTEMPTS = 1
"""int: How many times a re-roll is drawn. The generator prunes every wide binary to the rule (GEN.228), so one draw is enough; a system that still fails after it is reported, not looped on."""


def _problems(system):
    return [f"{problem.body}: {problem.message}" for problem in validation.check_star_system(system)]


REPORT_PROBLEMS = 3
"""int: How many of a system's broken bodies the console line names; the rest are counted."""


def _broke(item):
    """The bodies of a system that broke the rule, as the console line names them: the first few and the count of the rest."""
    lines = item["problems"]
    if not lines:
        return "unreadable"
    more = f"; and {len(lines) - REPORT_PROBLEMS} more" if len(lines) > REPORT_PROBLEMS else ""
    return "; ".join(lines[:REPORT_PROBLEMS]) + more


def _label(row):
    return f"{row['name']!r} (id {row['id']})"


def unstable(conn, say=None):
    """
    The wide binaries that fail the stability rule.

    Args:
        conn (store.Connection): An open connection to the galaxy.
        say (callable, optional): Takes one line; told how far the scan is.

    Returns:
        dict: `scanned` (the wide binaries read) and `found`, a list of dicts in row order: `id`, `name`, `sector_id`, `problems` (lines: the bodies that break the rule and how),
        `pair` (the pair's separation and closest approach, AU, from the stored row, or `None` where it could not be read), and `unreadable` (why the system could not be loaded, else `None`).
    """
    rows = conn.execute("SELECT id, name, sector_id FROM star_systems WHERE binary_configuration = 'wide' ORDER BY id").fetchall()
    found = []
    for number, row in enumerate(rows, 1):
        try:
            system = store.load_star_system(conn, row["id"])
            problems, pair, unreadable = _problems(system), _pair(system), None
        except Exception as exc:  # noqa: BLE001 -- a system that will not load cannot be told stable
            problems, pair, unreadable = [], None, f"{type(exc).__name__}: {exc}"
        if problems or unreadable:
            found.append({"id": row["id"], "name": row["name"], "sector_id": row["sector_id"], "problems": problems, "pair": pair,
                          "unreadable": unreadable})
        if say is not None and number % PROGRESS_EVERY == 0:
            say(f"Wide binaries: {number:,} of {len(rows):,} read, {len(found):,} break the rule so far.")
    return {"scanned": len(rows), "found": found}


def _pair(system):
    """`{"separation_au", "periapsis_au"}` of the pair of a loaded wide binary."""
    pair = getattr(system, "wide_binary", None)
    if pair is None:
        return None
    return {"separation_au": pair.separation_au, "periapsis_au": pair.periapsis_au}


def describe(system):
    """What a system made by the re-roll came out as, one phrase: the two stars' masses, the pair's closest approach and how many planets each star has."""
    primary, secondary = system.primary_star.mass / constants.SOLAR_MASS_TO_KG, system.secondary_star.mass / constants.SOLAR_MASS_TO_KG
    pair = system.wide_binary
    return (f"stars of {primary:.3g} and {secondary:.3g} solar masses, {pair.separation_au:.4g} AU apart (closest {pair.periapsis_au:.4g} AU), "
            f"{len(system.planets)} and {len(system.secondary_planets)} planets")


def _seed_address(conn, system_id):
    row = conn.execute("SELECT uid FROM star_systems WHERE id = ?", (system_id,)).fetchone()
    raw = None if row is None else row["uid"]
    return object_uid.format_id(object_uid.from_bytes(bytes(raw))) if raw is not None else f"row-{system_id}"


def blockers(conn, item, edited_sector_ids=()):
    """Why a system must be left alone, one phrase, or `None`: it hosts facilities, is built around a pre-placed bright star, or sits in a sector an admin edited."""
    found = store.system_content_blockers(conn, item["id"])
    if found["bright_star"]:
        return "it is built around a bright star the galaxy scatter placed"
    if found["facilities"]:
        return f"it hosts {found['facilities']} facilit{'y' if found['facilities'] == 1 else 'ies'}, which a re-roll would delete"
    if item["sector_id"] in set(edited_sector_ids):
        return "its sector was edited by an admin"
    return None


def reroll(conn, item):
    """
    Makes one unstable wide binary again, from its own seed (see the module's text), inside the caller's transaction.

    Args:
        conn (store.Connection): An open write connection.
        item (dict): One of `unstable`'s `found`.

    Returns:
        dict: `redone` (bool), `made` (`describe` of the new system, or `None`), `problems` (what the new system still breaks; empty when it follows the rule).
    """
    seed = store.get_galaxy_seed(conn) or mass_chain.FALLBACK_SEED
    address = _seed_address(conn, item["id"])
    skeleton = store.get_galaxy_shape(conn)
    shape_args = {} if skeleton is None else {"disk_scale_length_pc": skeleton.shape.disk_scale_length_pc, "shape": skeleton.shape}
    system = problems = None
    for _attempt in range(MAX_ATTEMPTS):
        config = SystemConfig()
        config.BINARY_SYSTEM = True
        config.WIDE_BINARY = True
        dist_ly = sector_context.apply_sector_context(conn, item["sector_id"], config) if item["sector_id"] is not None else None
        with galaxy_seed.seeded(seed, SEED_KIND, address, **shape_args):
            system = StarSystem(system_config=config, galactic_center_dist_ly=dist_ly)
            problems = _problems(system) if system.binary_type == "wide" else ["it did not come out a wide binary"]
        if not problems:
            break
    if problems:
        return {"redone": False, "made": describe(system) if system.binary_type == "wide" else None, "problems": problems}
    store.replace_star_system_content(conn, item["id"], system, config)
    return {"redone": True, "made": describe(system), "problems": []}


def run(config, edited_sector_ids=(), dry_run=False, say=print):
    """
    Finds the unstable wide binaries and re-rolls each, saying what it does (the console report Boss asked for: the system, the bodies that broke the rule, and what the re-roll made).

    Args:
        config (store.MySQLConfig): The galaxy database.
        edited_sector_ids (iterable): Sectors an admin edited, whose systems are left alone.
        dry_run (bool): Say which would be redone and change nothing.
        say (callable): Takes one line.

    Returns:
        dict: `scanned`, `unstable`, `redone`, `left` (systems left alone, with the reason) and `failed` (systems whose re-roll raised or came out unstable again, with the reason).
    """
    result = {"scanned": 0, "unstable": 0, "redone": 0, "left": [], "failed": []}
    conn = store.get_connection(config)
    try:
        scan = unstable(conn, say)
        result["scanned"], result["unstable"] = scan["scanned"], len(scan["found"])
        say(f"Wide binaries: {scan['scanned']:,} read, {len(scan['found']):,} break the stability rule.")
        for item in scan["found"]:
            why = f"could not be read ({item['unreadable']})" if item["unreadable"] else None
            why = why or blockers(conn, item, edited_sector_ids)
            broke = _broke(item)
            if why:
                result["left"].append((_label(item), why))
                say(f"star system {_label(item)}: broke the rule ({broke}); left alone: {why}.")
                continue
            if dry_run:
                say(f"star system {_label(item)}: broke the rule ({broke}); would be re-rolled.")
                continue
            try:
                with conn:
                    done = reroll(conn, item)
            except Exception as exc:  # noqa: BLE001 -- one system failing must not stop the others
                result["failed"].append((_label(item), str(exc)))
                say(f"star system {_label(item)}: broke the rule ({broke}); the re-roll FAILED: {exc}")
                continue
            if done["redone"]:
                result["redone"] += 1
                say(f"star system {_label(item)}: broke the rule ({broke}); re-rolled, now {done['made']}.")
            else:
                result["failed"].append((_label(item), "; ".join(done["problems"])))
                say(f"star system {_label(item)}: broke the rule ({broke}); the re-roll still breaks it ({'; '.join(done['problems'])}), left as it was.")
        say(f"Wide binaries: {result['redone']} re-rolled, {len(result['left'])} left alone, {len(result['failed'])} failed.")
        return result
    finally:
        conn.close()
