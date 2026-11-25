"""
planetgen/api/lookup.py

`GET /api/lookup?text=` (GEN.1140): what a printed Object ID or a name is, now. An ID or a name an object used to have answers with the object
that has it now (`docs/design/architecture/designation-history.md` section 5); a retired ID of a deleted object answers 410. The rules are
`galaxy/designations.resolve_uid` and `resolve_name`; this module reads the stored holders and the retired rows for them.
"""
import json
import time

from planetgen.db import designation_rows, local_name_search, store
from planetgen.galaxy import designations, object_uid, uid as sector_uids
from planetgen.names import naming_key
from planetgen.queue import api_jobs, redisqueue
from planetgen.util.format import format_duration_seconds
from . import ids

KIND_OF_TABLE = {table: kind for kind, table in ids.OBJECT_TABLES.items()}
KIND_OF_TABLE["sectors"] = ids.SECTOR
"""dict: The kind (the word before the colon of a reference) of each table that holds named objects."""

HOLDERS_PER_TABLE = 3
"""int: The most current holders of one name read from one table; a name is not unique among bodies."""

ALSO_LIMIT = 10
"""int: The most other past holders of a name listed in `also`."""

QUICK_SECONDS = 2.0
"""float: How long a bare local name is searched inside the request (API.38): the sectors nearest the person's own first. A search that is not done by then goes on as a job."""

PROGRESS_SECONDS = 1.0
"""float: The least time between two writes of a name search job's progress."""

PROGRESS_KEEP_SECONDS = 3600
"""int: How long a name search job's progress record is kept after its last write."""

SECONDS_PER_SECTOR = 0.00005
"""float: The time a job is expected to take for each sector it searches, before it has measured its own pace (`expected_seconds` of the job). Measured on the build box: the scramble is
0.009 ms a sector and the reads of a batch cost about as much again; this rounds up so the job is on the right lane."""


def _holders_of_name(conn, name):
    """`[(table, row id)]` of the objects named exactly `name` now, sectors and systems first."""
    holders = []
    for table in designation_rows.NAMED_TABLES:
        rows = conn.execute(f"SELECT id FROM {table} WHERE name = ? ORDER BY id LIMIT {HOLDERS_PER_TABLE}", (name,)).fetchall()
        holders.extend((table, row["id"]) for row in rows)
    return holders


def _object_answer(conn, table, row_id):
    """`{kind, ref, id, name, former_designations}` of one existing object, or `None` when the row is gone."""
    row = conn.execute(f"SELECT name FROM {table} WHERE id = ?", (row_id,)).fetchone()
    if row is None:
        return None
    kind = KIND_OF_TABLE[table]
    printed = ids.printed(conn, kind, row_id)
    return {"kind": kind, "ref": f"{kind}:{printed}", "id": printed, "name": row["name"],
            "former_designations": designation_rows.former_designations(conn, table, row_id)}


def lookup(conn, text, near=None, config=None):
    """
    Resolves a printed ID or a name.

    A bare local name (GEN.1255, `Zaidonr Jaisonn`, typed with no sector) that no object has as its stored name stands for one Object ID in every sector; it is searched sector by sector,
    the sectors nearest `near` first (`local_name_search`, API.38). The answer is at once when that takes `QUICK_SECONDS`, else 202 with the matches of the nearest sectors so far and a job
    (`GET /api/lookup/jobs/<id>`) that searches every sector and counts them as it goes.

    Returns:
        tuple: `(http status, body)`. 200 with `status` `current` (the object has the ID or name now) or `redirect` (it used to, and the body holds
            the object that has it now, with `retired` naming what was asked for); 410 `gone` for a designation of a deleted object; 404
            `unknown` when nobody ever had it. `also` lists other objects that used to have the same name, newest first.

    Args:
        near (str, optional): The 14 hex digits of the sector the person is in (a sector's designation), where a bare local name's search starts.
        config (MySQLConfig, optional): The database `conn` reads, for the job that carries a long search on.

    Raises:
        ValueError: Empty text, or `near` that is not a sector designation.
        planetgen.queue.api_jobs.NoQueue: A search that has to go on as a job and no queue to run it.
    """
    text = str(text or "").strip()
    if not text:
        raise ValueError("lookup needs text: an Object ID or a name")
    try:
        uid = object_uid.parse_id(text.upper())
    except ValueError:
        uid = None
    if uid is not None:
        return _by_id(conn, text, uid)
    typed = naming_key.uid_of_typed_name(text)     # GEN.1256: `Banhunm Finkap, sector 00FA0FFFE9E240`
    if typed is not None:
        status, body = _by_id(conn, text, typed)
        if status != 404:
            return status, body
    status, body = _by_name(conn, text)
    if status != 404:
        return status, body
    found = _by_local_name(conn, text, near, config)
    return (status, body) if found is None else found


def _by_id(conn, text, uid):
    retired = designation_rows.retired_uid(conn, uid)
    found = designations.resolve_uid(uid, lambda wanted: designation_rows.live_holder(conn, wanted), [] if retired is None else [retired])
    return _answer(conn, text, found, [])


def _by_name(conn, text):
    retired = designation_rows.retired_names(conn, text)
    name = retired[0].designation if retired else text
    found = designations.resolve_name(name, _holders_of_name(conn, name), retired)
    return _answer(conn, text, found, retired)


def _answer(conn, text, found, retired):
    if found.status == designations.UNKNOWN:
        return 404, {"query": text, "status": designations.UNKNOWN, "error": f"nothing has, or ever had, the ID or name {text}"}
    if found.status == designations.GONE:
        return designations.STATUS_GONE, {"query": text, "status": designations.GONE,
                                          "error": f"the object that had {text} was deleted", "gone_epoch": found.formerly[0].gone_epoch}
    body = _object_answer(conn, found.object_kind, found.object_row_id)
    if body is None:     # the row went without its tombstone: the same answer as a deleted object
        return designations.STATUS_GONE, {"query": text, "status": designations.GONE, "error": f"the object that had {text} was deleted"}
    body = {"query": text, "status": found.status, **body}
    if found.status == designations.REDIRECT:
        body["retired"] = text
    also = []
    for row in retired:
        if (row.object_kind, row.object_row_id) != (found.object_kind, found.object_row_id) and row.gone_epoch is None and len(also) < ALSO_LIMIT:
            other = _object_answer(conn, row.object_kind, row.object_row_id)
            if other is not None and other["ref"] != body["ref"]:
                also.append({"kind": other["kind"], "ref": other["ref"], "name": other["name"], "retired_epoch": row.retired_epoch})
    if also:
        body["also"] = also
    return 200, body


# -- a bare local name (API.38) --------------------------------------------------------------------------------------------------

def near_of(near):
    """The designation (an integer) of `near`, 14 hex digits, or `None` for no sector.

    Raises:
        ValueError: Text that is not a sector's designation.
    """
    if near is None or not str(near).strip():
        return None
    try:
        designation = int(str(near).strip().upper().removeprefix("SECTOR").strip(), 16)
        sector_uids.sector_address(designation)
    except ValueError as exc:
        raise ValueError(f"near is the 14 hex digits of a sector's designation, not {near!r}") from exc
    return designation


def _sectors_of(conn, near):
    """`(sector designations nearest first, designation of near or None)` for a pass over every sector that can hold an object named by its ID."""
    info = store.get_galaxy_shape(conn, light=True)
    edge_pc = getattr(info, "edge_pc", None) or 1.0        # the centres scale with the edge, so the order is the same for any edge
    return local_name_search.ordered(local_name_search.cells(conn), near_of(near), edge_pc)


def _by_local_name(conn, text, near, config):
    """
    `(status, body)` for a typed name that reads as a local name (`local_name_search.parts_of`), or `None` for any other text or when no naming key is in force.
    """
    key = naming_key.active_key()
    parts = local_name_search.parts_of(text, key)
    if parts is None:
        return None
    designations = _sectors_of(conn, near)
    started = time.monotonic()
    done = local_name_search.search(conn, parts, key, designations, deadline=started + QUICK_SECONDS)
    if done["complete"]:
        return _search_answer(conn, text, done)
    expected = max(1.0, (done["total"] - done["searched"]) * SECONDS_PER_SECTOR)
    job_id = api_jobs.submit(name_search_job, text, near, key, config, expected_seconds=expected)
    _write_progress(job_id, searched=0, total=done["total"], started=time.time(), updated=time.time(), partial=True)
    body = _search_body(conn, text, done)
    body.update(status="searching", job_id=job_id, status_url=f"/api/lookup/jobs/{job_id}",
                line=_line(done["searched"], done["total"], time.monotonic() - started))
    return 202, body


def _matches_of(conn, found):
    """`[{kind, ref, id, name, former_designations, sector}]` of the objects `local_name_search.search` found, in the order of the search; a row that went meanwhile is left out."""
    out = []
    for match in found["matches"]:
        answer = _object_answer(conn, match["table"], match["id"])
        if answer is not None:
            answer["sector"] = match["sector"]
            out.append(answer)
    return out


def _search_body(conn, text, found):
    matches = _matches_of(conn, found)
    return {"query": text, "status": "matches", "matches": matches, "searched": found["searched"], "of": found["total"], "complete": found["complete"]}


def _search_answer(conn, text, found):
    """`(status, body)` of a finished search: one match opens as `current` would, several are a list, none is `unknown` that says how many sectors were searched."""
    body = _search_body(conn, text, found)
    matches = body["matches"]
    if not matches:
        return 404, {"query": text, "status": designations.UNKNOWN, "searched": found["searched"], "of": found["total"], "complete": True,
                     "error": f"nothing has, or ever had, the ID or name {text}; every one of the {found['total']:,} sectors that can hold an object named by its ID was searched"}
    if len(matches) == 1:
        single = {"query": text, "status": designations.CURRENT, **matches[0], "searched": found["searched"], "of": found["total"], "complete": True}
        return 200, single
    return 200, body


def _line(searched, total, seconds):
    """The progress line of a search, with real counts and a clock time: `Searched 120,000 of 480,000 sectors in 1m 4s`."""
    return f"Searched {searched:,} of {total:,} sectors in {format_duration_seconds(seconds)}"


def _progress_key(job_id):
    return f"planetgen:namesearch:{job_id}"


def _write_progress(job_id, url=None, **fields):
    connection = redisqueue.connect(url or redisqueue.redis_url())
    key = _progress_key(job_id)
    connection.hset(key, mapping={name: json.dumps(value) for name, value in fields.items()})
    connection.expire(key, PROGRESS_KEEP_SECONDS)


def name_search_job(text, near, key, config):
    """
    The queued body of a bare local name's search that did not finish inside the request (API.38; a module-level function so a worker can import it by name): searches every sector
    nearest `near` first, writes its counts to the progress record (`status`) as it goes and returns the finished answer.

    Returns:
        dict: The answer the request would have given (`status` `current`, `matches` or `unknown`) with `http_status`.
    """
    from planetgen.queue import lightjobs
    job_id = lightjobs.current_id()
    names = {"key": key}
    naming_key.set_resolver(lambda: names["key"])
    conn = store.get_connection(config or store.DEFAULT_MYSQL_CONFIG)
    try:
        parts = local_name_search.parts_of(text, key)
        designations = _sectors_of(conn, near)
        started = time.monotonic()
        last = [0.0]

        def progress(searched, total):
            now = time.monotonic()
            if searched < total and now - last[0] < PROGRESS_SECONDS:
                return
            last[0] = now
            elapsed = now - started
            left = None if not searched else elapsed / searched * (total - searched)
            _write_progress(job_id, searched=searched, total=total, elapsed_seconds=elapsed, eta_seconds=left, line=_line(searched, total, elapsed), updated=time.time())

        found = local_name_search.search(conn, parts, key, designations, progress=progress)
        status, body = _search_answer(conn, text, found)
        body["http_status"] = status
        return body
    finally:
        conn.close()


def search_status(job_id):
    """
    Where a name search job stands, for `GET /api/lookup/jobs/<id>` (public: it names only the name searches): `state` (`queued`, `running`, `succeeded`, `failed`), `searched` of
    `of` sectors, `line` (the counts and the time in clock form), `eta_seconds` and `eta` (`1h 23m 45s`), and the finished answer as `result` once it ends. `None` for an id that is not a
    name search or has expired.

    Raises:
        planetgen.queue.api_jobs.NoQueue: No Redis server answers.
    """
    if not isinstance(job_id, str) or not api_jobs.lightjobs.JOB_ID_RE.match(job_id):
        return None
    try:
        record = redisqueue.connect(redisqueue.redis_url()).hgetall(_progress_key(job_id))
    except Exception as exc:  # noqa: BLE001 -- no Redis: the caller answers 503
        raise api_jobs.NoQueue(f"no Redis server answers: {exc}") from exc
    if not record:
        return None
    seen = {(name.decode() if isinstance(name, bytes) else name): json.loads(value) for name, value in record.items()}
    job = api_jobs.status(job_id) or {"state": "failed", "result": None, "error": "the job is gone"}
    body = {"id": job_id, "state": job["state"], "searched": seen.get("searched", 0), "of": seen.get("total"), "line": seen.get("line"),
            "eta_seconds": seen.get("eta_seconds"), "eta": None if seen.get("eta_seconds") is None else format_duration_seconds(seen["eta_seconds"]),
            "result": job["result"], "error": job["error"]}
    if job["state"] == "succeeded" and job["result"] is not None:
        body["searched"], body["of"] = job["result"].get("searched", body["searched"]), job["result"].get("of", body["of"])
    return body
