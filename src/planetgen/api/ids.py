# planetgen/api/ids.py

"""
Object IDs as the API's public reference (API.23).

Row ids (`star_systems.id` and the rest) stay the foreign and join keys inside the database and the
code that queries it; what the API, the pages, the URLs and `objectref` carry is the object's ID:
the 80-bit `uid` column printed as `FE81000A2B-0000005-000` (`galaxy/object_uid.py`), or for a sector
its designation printed as hex (`galaxy/uid.py`). This module is the one place that turns one into
the other:

- `row_id` / `row_ids` look a printed ID up in the object tables, `printed` / `printed_many` go the
  other way, in batches;
- `UidConverter` (the `uid` URL converter) accepts the printed forms in a route and returns the text
  upper-cased, so a page route needs no database to carry one;
- `resolve_view_args` (the app's `before_request` hook for API routes) turns the printed IDs in a
  route's URL into row ids before the view runs, so every view keeps working with row ids;
- `translate` rewrites the row ids in a JSON answer into printed IDs, at the places `ID_PATHS` lists
  for the endpoint, in one query per kind.

A kind prefix (`planet:FE81...`) is only a hint: the IDs of the object tables are unique on their own.
"""

import re
from collections import defaultdict

from werkzeug.routing import BaseConverter, ValidationError

from planetgen.galaxy import object_uid, uid as sector_uid

SECTOR = "sector"
"""str: The kind of a sector, whose ID is a plain integer (its designation)."""

OBJECT_TABLES = {
    "system": "star_systems",
    "star": "stars",
    "planet": "planets",
    "moon": "moons",
    "belt": "asteroid_belts",
    "comet": "comets",
    "facility": "facilities",
    "nebula": "nebulae",
    "asteroid_field": "asteroid_fields",
    "black_hole": "black_holes",
    "neutron_star": "neutron_stars",
    "supernova_remnant": "supernova_remnants",
    "rogue_planet": "rogue_planets",
    "interstellar_comet": "interstellar_comets",
    "quasar": "quasars",
}
"""dict: Each kind with an 80-bit `uid` column and the table that holds it. A phenomenon's kind is its type."""

PHENOMENON_KINDS = ("nebula", "asteroid_field", "black_hole", "neutron_star", "supernova_remnant",
                    "rogue_planet", "interstellar_comet", "quasar")
"""tuple: The kinds a generic `phenomenon` rule can resolve to (the row's `type`)."""

PRINTED_PATTERN = r"[0-9A-Fa-f]+(?:-[0-9A-Fa-f]+){0,2}"
"""str: What a printed ID looks like in a URL (the exact digit counts are checked by `parse`)."""

_SECTOR_RE = re.compile(r"^[0-9A-F]{1,16}$")
CHUNK = 1000
"""int: Ids looked up per query."""


class IdError(ValueError):
    """Text that is not a printed ID."""


def _table(kind):
    if kind == SECTOR:
        return "sectors"
    try:
        return OBJECT_TABLES[kind]
    except KeyError:
        raise IdError(f"unknown object kind: {kind!r}") from None


def parse(kind, text):
    """
    The stored value a printed ID names: an int for a sector, the 10 bytes for any other kind.

    Raises:
        IdError: For text that is not that kind's printed ID.
    """
    cleaned = str(text).strip().upper()
    if kind == SECTOR:
        if not _SECTOR_RE.match(cleaned):
            raise IdError(f"a sector ID is hex digits, got {text!r}")
        return int(cleaned, 16)
    _table(kind)
    try:
        return object_uid.to_bytes(object_uid.parse_id(cleaned))
    except ValueError as exc:
        raise IdError(str(exc)) from None


def format_stored(kind, value):
    """The printed form of a stored `uid` value (an int for a sector, bytes otherwise); `None` stays `None`."""
    if value is None:
        return None
    if kind == SECTOR:
        return sector_uid.format_sector_uid(int(value))
    return object_uid.format_id(object_uid.from_bytes(bytes(value)))


def looks_printed(text):
    """Whether `text` could be a printed ID (any kind)."""
    return bool(re.fullmatch(PRINTED_PATTERN, str(text).strip()))


def row_ids(conn, kind, printed_ids):
    """`{printed: row id}` for the printed IDs of `kind` that exist (the spelling is normalised to the
    upper-case printed form). Text that is not an ID is left out like one that names nothing."""
    table = _table(kind)
    wanted = {}
    for text in printed_ids:
        try:
            wanted[parse(kind, text)] = format_stored(kind, parse(kind, text))
        except IdError:
            continue
    found = {}
    values = list(wanted)
    for start in range(0, len(values), CHUNK):
        chunk = values[start:start + CHUNK]
        marks = ", ".join("?" for _ in chunk)
        for row in conn.execute(f"SELECT id, uid FROM {table} WHERE uid IN ({marks})", chunk).fetchall():
            stored = int(row["uid"]) if kind == SECTOR else bytes(row["uid"])
            found[wanted[stored]] = row["id"]
    return found


def row_id(conn, kind, printed_id):
    """The row id of the object a printed ID names, or `None`."""
    try:
        text = format_stored(kind, parse(kind, printed_id))
    except IdError:
        return None
    return row_ids(conn, kind, [text]).get(text)


def printed_many(conn, kind, ids):
    """`{row id: printed ID}` for the rows of `kind` among `ids`. A sector with no stored ID gets the ID it would
    be given (its designation from its address, else `uid.unplaced_sector_uid`)."""
    table = _table(kind)
    unique = sorted({int(i) for i in ids if i is not None})
    out = {}
    for start in range(0, len(unique), CHUNK):
        chunk = unique[start:start + CHUNK]
        marks = ", ".join("?" for _ in chunk)
        if kind == SECTOR:
            rows = conn.execute(
                f"SELECT id, uid, ring_index, layer_index, ring_slot_index FROM sectors WHERE id IN ({marks})",
                chunk).fetchall()
            for row in rows:
                out[row["id"]] = format_stored(SECTOR, _sector_value(row))
        else:
            for row in conn.execute(f"SELECT id, uid FROM {table} WHERE id IN ({marks})", chunk).fetchall():
                out[row["id"]] = format_stored(kind, row["uid"])
    return out


def _sector_value(row):
    if row["uid"] is not None:
        return row["uid"]
    if row["ring_index"] is not None:
        return sector_uid.sector_uid(row["ring_index"], row["layer_index"], row["ring_slot_index"])
    return sector_uid.unplaced_sector_uid(row["id"])


def printed(conn, kind, row_id_):
    """The printed ID of one row, or `None` when there is no such row."""
    return printed_many(conn, kind, [row_id_]).get(int(row_id_))


_REF_FORM = re.compile(r"^(?:([a-z_]+):)?(" + PRINTED_PATTERN + ")$")


def resolve_ref(conn, text):
    """
    The row-id reference (`system:12`, `objectref.parse`'s form) of a public reference: `<kind>:<printed ID>`
    (a bare three-part ID is a system's). A kind the API does not speak IDs for yet still takes its row id.

    Returns:
        str or None: `None` when the reference is well formed but names nothing.

    Raises:
        IdError: For text that is not a reference.
    """
    from planetgen.galaxy import objectref

    match = _REF_FORM.match(str(text).strip())
    if not match:
        raise IdError(f"not an object reference: {text!r}")
    kind, written = match.group(1), match.group(2).upper()
    if kind is None:
        if "-" not in written:
            raise IdError(f"not an object reference: {text!r}")
        kind = "system"
    if kind not in objectref.KINDS:
        raise IdError(f"unknown object kind: {kind!r}")
    if kind in ACTIVE_KINDS:
        parse(kind, written)  # IdError for text that cannot be this kind's ID
        found = row_id(conn, kind, written)
    elif written.isdigit():
        found = int(written)
    else:
        raise IdError(f"a {kind} is written with its row number for now, got {text!r}")
    return None if found is None else objectref.format(kind, found)


class UidConverter(BaseConverter):
    """The `uid` URL converter: a printed ID of any kind, upper-cased. It checks the shape only, so a page route
    needs no database; the API resolves it with `resolve_view_args`."""

    regex = PRINTED_PATTERN

    def to_python(self, value):
        return value.upper()

    def to_url(self, value):
        text = str(value).strip().upper()
        if not looks_printed(text):
            raise ValidationError(f"not an object ID: {value!r}")
        return text


# URL arguments of the API that carry an ID, and the kind each names. A `phenomenon_id` takes the kind of
# the route's `phenomenon_type` argument.
VIEW_ARG_KINDS = {
    "sector_id": SECTOR,
    "system_id": "system",
    "phenomenon_id": "phenomenon",
    "nebula_id": "nebula",
}


def resolve_view_args(conn_factory, view_args, abort_not_found, fallback=None):
    """
    Replaces the printed IDs in `view_args` (in place) with row ids. `conn_factory()` opens the request's
    read-only connection when an ID has to be looked up; `abort_not_found(message)` is called, and must
    raise, for an ID that names nothing; when it returns, `fallback` stands in.
    """
    for name, kind in VIEW_ARG_KINDS.items():
        if name not in view_args:
            continue
        if kind == "phenomenon":
            kind = view_args.get("phenomenon_type")
        text = view_args[name]
        found = row_id(conn_factory(), kind, text) if kind in ACTIVE_KINDS else None
        if found is None:
            abort_not_found(f"no {kind or 'object'} with ID {text}")
            found = fallback
        view_args[name] = found


ACTIVE_KINDS = frozenset(("sector", "system", *PHENOMENON_KINDS))
"""frozenset: The kinds whose IDs the API already speaks; any other kind keeps its row id for now."""

KEY_KINDS = {
    "sector_id": "sector",
    "star_system_id": "system",
    "system_id": "system",
    "capital_system_id": "system",
    "ref": "ref",
}
"""dict: Keys that name the same kind wherever they appear in an answer that has rules."""

_REF_TEXT = re.compile(r"^([a-z_]+):(\d+)$")
_NAV_KEY = re.compile(r"^phenomenon:([a-z_]+):(\d+)$")


def _kind_of(kind, container, context):
    """A rule's kind. `phenomenon` is the row's own `type` (else the route's `phenomenon_type`) and `by-kind` is the
    row's own `kind`."""
    if kind == "phenomenon":
        if isinstance(container, dict) and isinstance(container.get("type"), str):
            return container["type"]
        return context.get("phenomenon_type")
    if kind == "by-kind":
        return container.get("kind") if isinstance(container, dict) else None
    return kind


def _targets(kind, value, container, context):
    """`[(kind, row id)]` the value at a rule names: nothing for a value this rule does not apply to."""
    if kind == "ref":
        match = _REF_TEXT.match(value) if isinstance(value, str) else None
        return [(match.group(1), int(match.group(2)))] if match else []
    if kind == "navkey":
        if isinstance(value, int) and not isinstance(value, bool):
            return [("system", value)]
        match = _NAV_KEY.match(value) if isinstance(value, str) else None
        return [(match.group(1), int(match.group(2)))] if match else []
    if isinstance(value, int) and not isinstance(value, bool):
        return [(_kind_of(kind, container, context), value)]
    return []


def _replacement(kind, value, container, maps, context):
    """The value with its row id replaced by the printed ID, else `value` itself."""
    found = _targets(kind, value, container, context)
    if not found:
        return value
    target_kind, row = found[0]
    text = maps.get(target_kind, {}).get(row)
    if text is None:
        return value
    if kind == "ref":
        return f"{target_kind}:{text}"
    if kind == "navkey" and isinstance(value, str):
        return f"phenomenon:{target_kind}:{text}"
    return text


def _rule(rules, here, key):
    return rules.get(here) or (KEY_KINDS.get(key) if isinstance(key, str) else None)


def _collect(node, path, rules, wanted, context):
    """Adds the `kind -> row ids` every rule match under `node` names to `wanted`."""
    if isinstance(node, dict):
        for key, value in node.items():
            here = f"{path}.{key}" if path else str(key)
            kind = _rule(rules, here, key)
            if kind is not None:
                for target_kind, row in _targets(kind, value, node, context):
                    if target_kind in ACTIVE_KINDS:
                        wanted[target_kind].add(row)
            keys_kind = rules.get(here + "{}")
            if keys_kind is not None and isinstance(value, dict):
                for raw in value:
                    for target_kind, row in _targets(keys_kind, int(raw) if str(raw).isdigit() else raw, value, context):
                        if target_kind in ACTIVE_KINDS:
                            wanted[target_kind].add(row)
            _collect(value, here, rules, wanted, context)
    elif isinstance(node, list):
        here = path + "[]"
        kind = rules.get(here)
        for value in node:
            if kind is not None:
                for target_kind, row in _targets(kind, value, node, context):
                    if target_kind in ACTIVE_KINDS:
                        wanted[target_kind].add(row)
            _collect(value, here, rules, wanted, context)


def _rebuild(node, path, rules, maps, context):
    """A copy of `node` with the rule matches replaced by printed IDs (containers are copied, leaves shared)."""
    if isinstance(node, dict):
        out = {}
        for key, value in node.items():
            here = f"{path}.{key}" if path else str(key)
            kind = _rule(rules, here, key)
            if kind is not None and not isinstance(value, (dict, list)):
                value = _replacement(kind, value, node, maps, context)
            else:
                value = _rebuild(value, here, rules, maps, context)
                keys_kind = rules.get(here + "{}")
                if keys_kind is not None and isinstance(value, dict):
                    value = {(_replacement(keys_kind, int(k) if str(k).isdigit() else k, value, maps, context)
                              if str(k).isdigit() or keys_kind == "navkey" else k): v for k, v in value.items()}
            out[key] = value
        return out
    if isinstance(node, list):
        here = path + "[]"
        kind = rules.get(here)
        return [_replacement(kind, value, node, maps, context)
                if kind is not None and not isinstance(value, (dict, list)) else _rebuild(value, here, rules, maps, context)
                for value in node]
    return node


def translate(data, rules, conn, context=None):
    """
    A copy of `data` (a JSON answer) with the row ids that `rules` name replaced by printed IDs, one query per kind.
    `conn` is an open connection, or a function that opens one when an ID has to be looked up. `rules` is `{path: kind}` where a path is dotted keys with `[]` for each list level (`items[].id`,
    `route.path[]`); a path ending in `{}` renames the keys of that object (they are ids too). A kind is an
    object kind, or `phenomenon` (the row's `type`, else the route's `phenomenon_type` from `context`),
    `by-kind` (the row's `kind`), `ref` (an object reference `kind:id`) or `navkey` (a system id or a
    `phenomenon:type:id` NAV node). The keys in `KEY_KINDS` are translated wherever they appear. Only the kinds in
    `ACTIVE_KINDS` are; an id that names no row is left as it is. `data` itself is not touched, so a cached
    answer stays in row ids.
    """
    context = context or {}
    wanted = defaultdict(set)
    _collect(data, "", rules, wanted, context)
    wanted = {kind: ids for kind, ids in wanted.items() if kind in ACTIVE_KINDS and ids}
    if not wanted:
        return data
    conn = conn() if callable(conn) else conn
    maps = {kind: printed_many(conn, kind, ids) for kind, ids in wanted.items()}
    return _rebuild(data, "", rules, maps, context)


def _api_request():
    from flask import request

    return request.path.startswith("/api/")


def _resolve(view_args, strict):
    """Replaces the printed IDs in `view_args` (in place) with row ids. With `strict` an ID that names nothing ends
    in a 404; otherwise row 0 (which names nothing) stands in."""
    from planetgen.api.common import ApiError
    from planetgen.api.routes import get_db

    def missing(message):
        if strict:
            raise ApiError(message, 404)

    resolve_view_args(get_db, view_args, missing, fallback=0)


def _resolve_request_ids():
    """The app's `before_request` hook (after the query-string and body checks): row ids for the printed IDs in an
    API route's URL. A route that needs an admin resolves its own, in `authz.require_admin` once the caller is
    known; a write resolves to row 0 when nothing has the ID, and the route answers 404 itself."""
    from flask import current_app, request

    view_args = request.view_args
    if not view_args or not _api_request():
        return
    view = current_app.view_functions.get(request.endpoint)
    if getattr(view, "admin_required", None) is not None:
        return
    _resolve(view_args, strict=request.method in ("GET", "HEAD"))


def resolve_kwargs(kwargs):
    """`require_admin`'s call, with the view's keyword arguments: their printed IDs become row ids (a 404 for one
    that names nothing)."""
    from flask import request

    if kwargs and _api_request():
        _resolve(kwargs, strict=True)


def translate_response(obj):
    """`obj` (a JSON answer about to be written) with the row ids of this endpoint's `ID_PATHS` printed; `obj`
    itself when the request is not an API one or its endpoint has no rules."""
    from flask import has_request_context, request

    from planetgen.api.idpaths import ID_PATHS

    if not has_request_context() or not _api_request():
        return obj
    rules = ID_PATHS.get(request.endpoint)
    if rules is None or not isinstance(obj, (dict, list)):
        return obj
    from planetgen.api.routes import get_db

    def connection():
        conn = get_db()
        if request.method not in ("GET", "HEAD"):
            conn.rollback()  # a write may have committed since this connection's snapshot began
        return conn

    context = {"phenomenon_type": (request.view_args or {}).get("phenomenon_type")}
    return translate(obj, rules, connection, context)


def install(app):
    """Registers the `uid` URL converter and the hook that resolves printed IDs in API routes. Call it after the
    app's own request checks, which must see a request before any database is asked."""
    app.url_map.converters["uid"] = UidConverter
    app.before_request(_resolve_request_ids)
