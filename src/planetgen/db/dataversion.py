# planetgen/db/dataversion.py

"""
Data versions for the shared page cache (PERF.104, docs/design/operations/cache-warm-up.md).

A cached answer names the pieces of data it read (its tags, `web/lib/pagetags.py`) and is
served only while every one of them still has the counter it had when the answer was made.
Whatever changes data bumps the counters of the tags that change touches (`record`), in
Redis, so every process sees it at once and a stale copy is never served as current.

One Redis hash per database, `planetgen:ver:<database>`: tag name to counter, and `_id`, an
incarnation token made when the hash is first written. A Redis that lost its data (a flush, a
cold start without the append-only file) starts a new incarnation, so a copy kept in a web
process from before it is never taken for current.

Tags (`TAG_*`, `sector_tag` and the like):

- `epoch`: every answer carries it; bumped by anything that cannot be placed more exactly.
- `sector:<id>`, `system:<id>`, `phenomenon:<type>:<id>`: one object's own page and what lists it (ids as printed, API.23).
- `sector-set`, `system-set`, `phenomena-set`: an object was added or deleted (who is next to whom, counts).
- `sector-data`, `system-data`, `phenomena-data`: any such object was added, deleted or edited (what a list shows).
- `system-any`: every system page (an edit to a star, planet, moon or facility, whose system the
  route does not name).
- `population`, `positions`: the population pass and the orbital update.

Nothing here raises for Redis: a writer that cannot reach it prints a warning (the readers
cannot reach it either, so they cache nothing); `versions` raises `Unavailable` and the page
cache then serves nothing from its copies.
"""

import secrets
import sys

import redis

from planetgen.queue import redisclient, warmstart

TAG_EPOCH = "epoch"
TAG_SECTOR_SET = "sector-set"
TAG_SYSTEM_SET = "system-set"
TAG_PHENOMENA_SET = "phenomena-set"
TAG_SECTOR_DATA = "sector-data"
TAG_SYSTEM_DATA = "system-data"
TAG_PHENOMENA_DATA = "phenomena-data"
TAG_SYSTEM_ANY = "system-any"
TAG_POPULATION = "population"
TAG_POSITIONS = "positions"

ID_FIELD = "_id"
"""str: The hash field that holds the incarnation token."""

COUNT_FIELD = "_n"
"""str: The hash field every bump adds one to, so a reader can tell that something changed while it fetched."""

STATE_FIELD = "_state"
"""str: The hash field that holds the last galaxy content state (`db.query.galaxy_changes`) whose changes were bumped."""

SET_TAGS = (TAG_SECTOR_SET, TAG_SYSTEM_SET, TAG_PHENOMENA_SET)
"""tuple: The tags that say an object was added or deleted."""

DATA_TAGS = (TAG_SECTOR_DATA, TAG_SYSTEM_DATA, TAG_PHENOMENA_DATA)
"""tuple: The tags that say an object of that kind was added, deleted or edited (what a list shows)."""


class Unavailable(Exception):
    """Redis does not answer, so no version can be read."""


def sector_tag(sector_id):
    """The tag of a sector, from its printed ID (API.23), the spelling the pages and the API use."""
    return f"sector:{str(sector_id).strip().upper()}"


def system_tag(system_id):
    return f"system:{str(system_id).strip().upper()}"


def phenomenon_tag(kind, phenomenon_id):
    return f"phenomenon:{kind}:{str(phenomenon_id).strip().upper()}"


def key(database):
    """The hash for `database`."""
    return f"planetgen:ver:{database}"


KEEP_SECONDS = 30 * 86400
"""int: A version hash nothing wrote to or read for this long goes (a database dropped, a test's); the answers
that used it are long gone by then (`page_cache.shared_seconds`)."""


def _client():
    return redisclient.client(timeout=0.5, retries=1, decode_responses=True)


def _ensure(client, database):
    """Makes the hash (with its incarnation token) when it is missing."""
    if client.hsetnx(key(database), ID_FIELD, secrets.token_hex(8)):
        client.expire(key(database), KEEP_SECONDS)


def incarnation(database):
    """The token of `database`'s version hash, made when missing. Raises `Unavailable`."""
    try:
        client = _client()
        _ensure(client, database)
        return client.hget(key(database), ID_FIELD)
    except (redis.RedisError, OSError) as exc:
        raise Unavailable(str(exc)) from exc


def versions(database, tags):
    """
    The counters of `tags` and the incarnation, in one round trip.

    Returns:
        tuple: `(incarnation, {tag: counter})`; a tag never bumped is 0.

    Raises:
        Unavailable: Redis does not answer.
    """
    tags = list(dict.fromkeys(tags))
    try:
        client = _client()
        _ensure(client, database)
        values = client.hmget(key(database), [ID_FIELD, *tags])
    except (redis.RedisError, OSError) as exc:
        raise Unavailable(str(exc)) from exc
    return values[0], {tag: int(value or 0) for tag, value in zip(tags, values[1:])}


def bump(database, *tags):
    """
    Adds one to each of `tags` for `database`, then asks for the warm-up (`warmstart.request`). Returns whether the
    change is recorded.

    A Redis that refuses the write (a full `noeviction` server answers OOM to anything that adds data) still
    allows a delete, so the whole version hash is dropped instead: a new incarnation, which makes every copy
    taken under the old one stale. That is a wider invalidation, never a stale page; it is warned about.
    """
    tags = list(dict.fromkeys(tags))
    if not tags or database is None:
        return True
    try:
        client = _client()
        pipe = client.pipeline(transaction=True)
        pipe.hsetnx(key(database), ID_FIELD, secrets.token_hex(8))
        pipe.hincrby(key(database), COUNT_FIELD, 1)
        for tag in tags:
            pipe.hincrby(key(database), tag, 1)
        pipe.expire(key(database), KEEP_SECONDS)
        pipe.execute()
    except (redis.RedisError, OSError) as exc:
        try:
            _client().delete(key(database))
        except (redis.RedisError, OSError):
            print(f"warning: the page cache's data versions were not updated ({exc}); "
                  "cached pages may be stale until Redis answers.", file=sys.stderr)
            return False
        print(f"warning: Redis refused the page cache's data versions ({exc}); every cached page was made stale "
              "instead. If Redis is out of memory, raise maxmemory.", file=sys.stderr)
    # Whatever the page cache cannot serve any more is rebuilt by the warm-up, which starts by itself (PERF.141).
    warmstart.request(database, tags[0])
    return True


def tags_for_change(kind, sector_id=None, system_id=None, phenomenon=None):
    """
    The tags a change of `kind` bumps (the table in docs/design/operations/cache-warm-up.md).

    Args:
        kind (str): One of `CHANGES`.
        sector_id (int, optional): The sector concerned.
        system_id (int, optional): The system concerned.
        phenomenon (tuple, optional): `(type, id)`.

    Raises:
        ValueError: An unknown kind, or an id the kind needs and was not given.
    """
    if kind not in CHANGES:
        raise ValueError(f"unknown kind of data change {kind!r}; one of {', '.join(sorted(CHANGES))}")
    tags = []
    if kind == "sector-edited":
        tags += [sector_tag(_need(sector_id, kind)), TAG_SECTOR_DATA]
    elif kind == "sector-contents":
        tags += [sector_tag(_need(sector_id, kind)), TAG_SECTOR_DATA, TAG_SYSTEM_SET, TAG_SYSTEM_DATA,
                 TAG_PHENOMENA_SET, TAG_PHENOMENA_DATA]
    elif kind == "sector-set":
        tags += [*SET_TAGS, *DATA_TAGS]
        if sector_id is not None:
            tags.append(sector_tag(sector_id))
    elif kind == "system-edited":
        tags += [system_tag(_need(system_id, kind)), TAG_SYSTEM_DATA]
        if sector_id is not None:
            tags.append(sector_tag(sector_id))
    elif kind == "system-any":
        tags += [TAG_SYSTEM_ANY, TAG_SYSTEM_DATA]
    elif kind == "system-set":
        tags += [TAG_SYSTEM_SET, TAG_SYSTEM_DATA, TAG_SECTOR_DATA]
        if system_id is not None:
            tags.append(system_tag(system_id))
        if sector_id is not None:
            tags.append(sector_tag(sector_id))
    elif kind == "phenomenon-edited":
        tags += [phenomenon_tag(*_need(phenomenon, kind)), TAG_PHENOMENA_DATA]
        if sector_id is not None:
            tags.append(sector_tag(sector_id))
    elif kind == "phenomena-set":
        tags += [TAG_PHENOMENA_SET, TAG_PHENOMENA_DATA]
        if phenomenon is not None:
            tags.append(phenomenon_tag(*phenomenon))
        if sector_id is not None:
            tags.append(sector_tag(sector_id))
    elif kind == "generated":
        tags += [*SET_TAGS, *DATA_TAGS]
    elif kind == "population":
        tags.append(TAG_POPULATION)
    elif kind == "orbits":
        tags.append(TAG_POSITIONS)
    else:  # everything
        tags.append(TAG_EPOCH)
    return tags


CHANGES = {
    "sector-edited": "a sector's own fields (a rename, a wiki link)",
    "sector-contents": "a sector's systems and phenomena were deleted or regenerated (its own and its systems' pages)",
    "sector-set": "a sector was added or deleted",
    "system-edited": "a system's own fields, or a star, planet, moon or facility in it",
    "system-any": "a star, planet, moon or facility edit whose system is not known",
    "system-set": "a system was added or deleted",
    "phenomenon-edited": "a phenomenon's own fields",
    "phenomena-set": "a phenomenon was added or deleted",
    "generated": "a generation run placed sectors, systems and phenomena",
    "population": "the population pass ran, or a species, polity or territory was edited",
    "orbits": "the orbital update moved objects",
    "everything": "the plan, the shape, the naming key, a release, or a write nothing above describes",
}
"""dict: Every kind of data change `record` knows, and what it is."""


def _need(value, kind):
    if value is None:
        raise ValueError(f"a {kind!r} change needs the id of what changed")
    return value


def record(database, kind, **ids):
    """
    Notes a change of `kind` in `database` (see `tags_for_change`): bumps its tags.

    Returns:
        list: The tags bumped.
    """
    tags = tags_for_change(kind, **ids)
    bump(database, *tags)
    return tags


def claim_state(database, state, previous):
    """
    Records that the galaxy's changes up to `state` (a token of `db.query.galaxy_changes`) have been bumped,
    when `previous` is still the recorded state. Several web processes see the same change; the first to claim
    it bumps, the others skip.

    Returns:
        bool: Whether this caller won and should bump. Raises `Unavailable` when Redis does not answer.
    """
    script = ("local cur = redis.call('HGET', KEYS[1], ARGV[1]); "
              "if (cur == false and ARGV[3] == '') or cur == ARGV[3] then "
              "redis.call('HSET', KEYS[1], ARGV[1], ARGV[2]); return 1 end; return 0")
    try:
        client = _client()
        _ensure(client, database)
        return bool(client.eval(script, 1, key(database), STATE_FIELD, state, previous or ""))
    except (redis.RedisError, OSError) as exc:
        raise Unavailable(str(exc)) from exc


def recorded_state(database):
    """The state `claim_state` recorded last, or `None`. Raises `Unavailable`."""
    try:
        return _client().hget(key(database), STATE_FIELD)
    except (redis.RedisError, OSError) as exc:
        raise Unavailable(str(exc)) from exc
