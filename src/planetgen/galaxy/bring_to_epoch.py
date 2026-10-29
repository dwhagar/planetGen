# planetgen/galaxy/bring_to_epoch.py

"""
Bringing a stored galaxy to the current generator epoch without a full regeneration (GEN.764, step 3).

`survey` reads what the galaxy's stages were made under (`epoch_record`) and says what a run would do; `run` does it:

- A galaxy whose plan is behind (an epoch since changed the galaxy's shape) cannot be brought up in part: the run
  refuses and says a full regeneration is needed, naming the epochs.
- Likewise a scatter that is behind: the sectors are filled from the scatters, so they are redrawn first, with
  `planetgen clear-scatter` and `planetgen scatter` (this command does not redraw a scatter yet).
- Each galaxy-placed sector whose fill is behind is generated again under the current code, in address order
  (`queue.api_jobs.regenerate_sector`: the sector is deleted and made again from its address and the galaxy's seed, then
  the paths around it are settled). The sector records its epoch when it is saved, so an interrupted run is resumed
  by running it again.
- A sector with admin edits is skipped and listed (Boss's choice, GEN.764): an edit has no row of its own, so a sector
  counts as edited when the admin audit log names an edit of the sector, one of its systems, bodies or facilities, or
  of a phenomenon filed in it. An edit of something since deleted cannot be traced to its sector.
- The species made under an older population epoch outside the sectors redone are made again (a species edited by hand,
  `tech_source = 'admin'`, is kept), then the population pass runs.
"""

import pymysql

from planetgen.db import kind_tables, store
from planetgen.galaxy import epoch_record, uid as sector_uid
from planetgen.population import model as population
from planetgen.util import generator_epoch, log

_CONTENT_KINDS = ("system", "star", "planet", "moon", "belt", "comet", "facility")
"""tuple: The reference kinds of an audit target (`system:12`) that sit in a system or a sector."""

_NOT_AN_EDIT = ("sector.regenerate", "sector.delete", "sector.create", "facility.suggest_ratings")
"""tuple: The audit actions that name a sector but do not leave an edit in it: a regeneration replaces what was there,
a deleted sector is gone, a created one is hand-made (no grid address) and the ratings are filled from the facility."""

PAGE_LINES = 20
"""int: How many sectors of each kind the Generate page names; the rest are counted."""

_BODY_TABLES = {"planet": "planets", "moon": "moons", "belt": "asteroid_belts", "comet": "comets", "star": "stars"}


def _sector_of(conn, kind, row_id):
    """The id of the sector the object `kind:row_id` is filed in, or `None` (no such row, or one in no sector)."""
    if kind == "sector":
        row = conn.execute("SELECT id FROM sectors WHERE id = ?", (row_id,)).fetchone()
        return None if row is None else row["id"]
    if kind == "system":
        row = conn.execute("SELECT sector_id FROM star_systems WHERE id = ?", (row_id,)).fetchone()
    elif kind in _BODY_TABLES:
        row = conn.execute(
            f"SELECT ss.sector_id FROM {_BODY_TABLES[kind]} b JOIN star_systems ss ON ss.id = b.star_system_id "
            "WHERE b.id = ?", (row_id,)).fetchone()
    elif kind == "facility":
        row = conn.execute("SELECT COALESCE(f.sector_id, ss.sector_id) AS sector_id FROM facilities f "
                           "LEFT JOIN star_systems ss ON ss.id = f.star_system_id WHERE f.id = ?", (row_id,)).fetchone()
    else:
        table = next((entry.table for entry in kind_tables.KIND_TABLES if entry.generator_type == kind), None)
        if table is None:
            return None
        row = conn.execute(f"SELECT sector_id FROM {table} WHERE id = ?", (row_id,)).fetchone()
    return None if row is None else row["sector_id"]


def edited_sectors(conn, control_conn):
    """
    The sectors the admin audit log says an admin edited (`authz.audit`: one row per edit, its target `kind:row id`).

    Returns:
        dict: `{sector row id: [audit action, ...]}`, oldest edit first; empty when there is no control connection.
    """
    if control_conn is None:
        return {}
    found = {}
    rows = control_conn.execute("SELECT action, target FROM admin_audit_log WHERE target IS NOT NULL ORDER BY id").fetchall()
    for row in rows:
        kind, _, text = row["target"].partition(":")
        if row["action"] in _NOT_AN_EDIT or not text.isdigit():
            continue
        sector_id = _sector_of(conn, kind, int(text))
        if sector_id is not None:
            found.setdefault(sector_id, []).append(row["action"])
    return found


def _label(row):
    return f"{sector_uid.format_sector_uid(sector_uid.sector_uid(row['shell_index'], row['band_index'], row['wedge_index']))} " \
           f"(shell {row['shell_index']} band {row['band_index']} wedge {row['wedge_index']})"


def _numbers(epochs):
    return ", ".join(str(epoch) for epoch in epochs)


def survey(conn, control_conn=None, current=None):
    """
    What bringing the galaxy to epoch `current` (default the running one) would do.

    Returns:
        dict: `current`; `behind` (`epoch_record.behind`); `refuse`, the reasons a run cannot go on (lines; empty
        when it can); `refill`, the galaxy-placed sectors to generate again (dicts of `id`, `label`, `epoch`), in
        address order; `skipped`, the sectors left alone (the same with `why`); `handmade`, how many stale sectors
        have no grid address; `species`, the stale species to make again (`redo`) and those kept as edited (`kept`).
    """
    current = generator_epoch.generator_epoch() if current is None else current
    found = epoch_record.behind(conn, current)
    refuse = []
    if "plan" in found["stale"]:
        epochs = generator_epoch.epochs_changing("plan", found["stale"]["plan"], current)
        refuse.append(f"The plan was made under epoch {found['stale']['plan']} and epoch {_numbers(epochs)} changed the "
                      "galaxy's shape: a full regeneration is needed, a part of the galaxy cannot be brought up.")
    for kind, epoch in found["stale"]["scatter"].items():
        epochs = generator_epoch.epochs_changing("scatter", 1 if epoch is None else epoch, current)
        refuse.append(f"The {kind} scatter was made under epoch {epoch} and epoch {_numbers(epochs)} changed it: "
                      f"draw it again (planetgen clear-scatter, then planetgen scatter) before bringing the sectors up.")
    edited = edited_sectors(conn, control_conn)
    refill, skipped, handmade = [], [], 0
    rows = conn.execute("SELECT id, shell_index, band_index, wedge_index, generator_epoch FROM sectors "
                        "ORDER BY shell_index, band_index, wedge_index, id").fetchall()
    for row in rows:
        if not found["stale"]["fill"] or not epoch_record.needs("fill", row["generator_epoch"], current):
            continue
        if row["shell_index"] is None:
            handmade += 1
            continue
        item = {"id": row["id"], "label": _label(row), "epoch": row["generator_epoch"]}
        if row["id"] in edited:
            actions = edited[row["id"]]
            skipped.append({**item, "why": f"edited by an admin ({', '.join(sorted(set(actions)))})"})
        else:
            refill.append(item)
    species = {"redo": 0, "kept": 0}
    if found["stale"]["population"]:
        where, values = _stale_species(found["stale"]["population"])
        counts = conn.execute(f"SELECT COALESCE(SUM(tech_source = 'admin'), 0) AS kept, COUNT(*) AS total "
                              f"FROM species WHERE {where}", values).fetchone()
        species = {"redo": int(counts["total"]) - int(counts["kept"]), "kept": int(counts["kept"])}
    return {"current": current, "behind": found, "refuse": refuse, "refill": refill, "skipped": skipped,
            "handmade": handmade, "species": species}


def status(config, control_config=None):
    """
    What `run` would do, for the Generate page (GEN.764, step 4), changing nothing.

    Args:
        config (MySQLConfig): The galaxy database.
        control_config (MySQLConfig, optional): The control database, whose audit log names the edited sectors.

    Returns:
        dict: `survey`'s dict with `lines` (`report_lines`), `nothing` (the galaxy is at the running epoch) and `error`
        (`None`), or only `error` (why the galaxy could not be read, or the audit log is unreadable so no edited sector
        could be told apart).
    """
    conn = store.get_connection(config, ensure_schema=False)
    control_conn = None
    try:
        if control_config is not None:
            try:
                control_conn = store.get_control_connection(control_config, ensure_schema=False)
            except pymysql.MySQLError as exc:
                return {"error": f"The admin audit log could not be read ({exc}); without it the edited sectors "
                                 "cannot be told apart."}
        found = survey(conn, control_conn)
        return {**found, "lines": report_lines(found, limit=PAGE_LINES), "nothing": _nothing_to_do(found), "error": None}
    finally:
        if control_conn is not None:
            control_conn.close()
        conn.close()


def _nothing_to_do(found):
    """Whether the galaxy has no record behind the epoch (so a run has nothing to say or to do)."""
    return not (found["refuse"] or found["refill"] or found["skipped"] or found["handmade"] or found["species"]["redo"]
                or any(found["behind"]["stale"].values()))


def report_lines(found, limit=None):
    """The lines that say what a run would do (or did not do), for a dry run and the start of a run. With `limit`, only
    that many sectors of each kind are named and the rest are counted (the Generate page)."""
    lines = []
    if _nothing_to_do(found):
        return [f"The galaxy is at epoch {found['current']}: nothing to bring up."]
    lines.append(f"Bringing the galaxy to epoch {found['current']}; stages behind: "
                 f"{', '.join(found['behind']['stages']) or 'none'}.")
    lines.extend(found["refuse"])
    if found["refuse"]:
        lines.append("Nothing was changed.")
        return lines
    for item in found["refill"][:limit]:
        lines.append(f"sector {item['label']}: made under epoch {item['epoch']}, to be generated again.")
    if limit is not None and len(found["refill"]) > limit:
        lines.append(f"... and {len(found['refill']) - limit:,} more sector(s) to be generated again.")
    for item in found["skipped"][:limit]:
        lines.append(f"sector {item['label']}: made under epoch {item['epoch']}, left alone: {item['why']}.")
    if limit is not None and len(found["skipped"]) > limit:
        lines.append(f"... and {len(found['skipped']) - limit:,} more sector(s) left alone.")
    if found["handmade"]:
        lines.append(f"{found['handmade']} hand-made sector(s) with no grid address are not generated and are left alone.")
    species = found["species"]
    if species["redo"] or species["kept"]:
        lines.append(f"{species['redo']} species made under an older epoch to be made again"
                     + (f", {species['kept']} edited by hand kept." if species["kept"] else "."))
    lines.append(f"{len(found['refill'])} sector(s) to redo, {len(found['skipped'])} left alone.")
    return lines


def _stale_species(epochs):
    """`(condition, values)` for the species made under one of `epochs` (a None is a species with no epoch)."""
    known = [epoch for epoch in epochs if epoch is not None]
    marks = ", ".join("?" for _ in known)
    parts = ([f"generator_epoch IN ({marks})"] if known else []) + (["generator_epoch IS NULL"] if None in epochs else [])
    return "(" + " OR ".join(parts) + ")", known


def _redo_species(conn, epochs):
    """Deletes the species made under one of the stale population `epochs` (not the hand-edited ones) and moves the
    population pass's watermark back to the lowest of their worlds, so the pass makes them again. Returns how many."""
    if not epochs:
        return 0
    where, values = _stale_species(epochs)
    where += " AND (tech_source IS NULL OR tech_source <> 'admin')"
    low = conn.execute(f"SELECT MIN(homeworld_planet_id) AS low FROM species WHERE {where}", values).fetchone()["low"]
    if low is None:
        return 0
    with conn:
        count = conn.execute(f"DELETE FROM species WHERE {where}", values).rowcount
        conn.execute("UPDATE population_state SET scanned_planet_id = LEAST(scanned_planet_id, ?) WHERE id = 1", (low - 1,))
    return count


def run(config, control_config=None, dry_run=False, say=log.normal):
    """
    Brings the stored galaxy to the running epoch (see the module's text), saying each step with `say`.

    Args:
        config (MySQLConfig): The galaxy database.
        control_config (MySQLConfig, optional): The control database, whose audit log names the edited sectors.
        dry_run (bool): Say what would be done and change nothing.
        say (callable): Takes one line.

    Returns:
        dict: `refused` (bool), `redone` (sectors generated again), `skipped` (sectors left alone), `failed` (sectors
        that raised, with the reason), `species` (made again).
    """
    from planetgen.queue import api_jobs
    result = {"refused": False, "redone": 0, "skipped": 0, "failed": [], "species": 0}
    conn = store.get_connection(config)
    control_conn = None
    if control_config is not None:
        try:
            control_conn = store.get_control_connection(control_config, ensure_schema=False)
        except pymysql.MySQLError as exc:
            conn.close()
            say(f"The admin audit log could not be read ({exc}); without it the edited sectors cannot be told apart "
                "and none is generated again. Nothing was changed.")
            result["refused"] = True
            return result
    try:
        found = survey(conn, control_conn)
        for line in report_lines(found):
            say(line)
        result["skipped"] = len(found["skipped"])
        if found["refuse"]:
            result["refused"] = True
            return result
        if dry_run or _nothing_to_do(found):
            return result
        for number, item in enumerate(found["refill"], 1):
            try:
                made = api_jobs.regenerate_sector(item["id"], config)
            except Exception as exc:  # noqa: BLE001 -- one sector failing must not stop the others
                result["failed"].append((item["label"], str(exc)))
                say(f"sector {item['label']} ({number} of {len(found['refill'])}): FAILED: {exc}")
                continue
            result["redone"] += 1
            say(f"sector {item['label']} ({number} of {len(found['refill'])}): generated again under epoch "
                f"{found['current']}, now {made['sector_name']}.")
        result["species"] = _redo_species(conn, list(found["behind"]["stale"]["population"]))
        if result["species"] or result["redone"]:
            say("Population: making the species again.")
            counts = population.run_population(conn)
            say(f"Population: {counts['new_species']:,} new species; {counts['species']:,} in all.")
            owned = population.run_territories(conn)
            say(f"Population: {owned['polities']:,} polities own {owned['owned_systems']:,} systems.")
        say(f"Done: {result['redone']} sector(s) generated again, {result['skipped']} left alone"
            + (f", {len(result['failed'])} failed (run it again to retry them)." if result["failed"] else "."))
        return result
    finally:
        if control_conn is not None:
            control_conn.close()
        conn.close()
