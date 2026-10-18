# planetgen/web/lib/pagetags.py

"""
What each cached answer reads, and what each write changes (PERF.104, docs/design/cache-warm-up.md).

`read_tags` names the data tags (`planetgen.db.dataversion`) an API answer depends on, from its
path and its body: a sector page reads its own sector, each neighbour it links to, the sets of
sectors and phenomena; a system page its own system and its sector; a list the sets it lists.
`change_for_write` maps a write through the API to the change kind (`dataversion.CHANGES`) it
makes. These two tables are the whole contract between the page cache and the writers; a path
neither knows reads everything and a write neither knows changes everything.
"""

import json
import re

from planetgen.db import dataversion as dv

_ID = r"([^/]+)"
_SYSTEM_SUBPATHS = ("text", "description", "scene", "facilities", "habitable-species", "owner")

ALL_TAGS = (*dv.SET_TAGS, *dv.DATA_TAGS, dv.TAG_SYSTEM_ANY, dv.TAG_POPULATION, dv.TAG_POSITIONS)
"""tuple: What an answer nothing here describes reads, besides `epoch`."""


def _body(body):
    try:
        return json.loads(body)
    except (TypeError, ValueError):
        return None


def _neighbour_tags(detail):
    tags = []
    for neighbour in (detail or {}).get("neighbors") or []:
        if isinstance(neighbour, dict) and neighbour.get("sector_id"):
            tags.append(dv.sector_tag(neighbour["sector_id"]))
    return tags


def read_tags(path, body):
    """
    The tags the answer `body` of `GET /api<path>` reads.

    Args:
        path (str): The path under `/api`, no query.
        body (str or bytes): The raw JSON body, read for the ids of the objects it links to.

    Returns:
        list: Tag names, `epoch` first.
    """
    tags = [dv.TAG_EPOCH]

    match = re.fullmatch(rf"/sectors/{_ID}(?:/facilities)?", path)
    if match:
        detail = _body(body) if path.count("/") == 2 else None
        return tags + [dv.sector_tag(match.group(1)), dv.TAG_SECTOR_SET, dv.TAG_PHENOMENA_DATA, dv.TAG_POPULATION,
                       dv.TAG_POSITIONS, *_neighbour_tags(detail)]

    match = re.fullmatch(rf"/systems/{_ID}(?:/({'|'.join(_SYSTEM_SUBPATHS)}))?", path)
    if match:
        detail = _body(body) if match.group(2) is None else None
        sector = (detail or {}).get("sector_id") if isinstance(detail, dict) else None
        tags += [dv.system_tag(match.group(1)), dv.TAG_SYSTEM_ANY, dv.TAG_PHENOMENA_DATA, dv.TAG_POPULATION, dv.TAG_POSITIONS]
        # The sector's tag when the body names it; otherwise every set, since the sector is not known.
        return tags + ([dv.sector_tag(sector)] if sector else [dv.TAG_SECTOR_SET, dv.TAG_SYSTEM_SET])

    match = re.fullmatch(rf"/phenomena/{_ID}/{_ID}", path)
    if match:
        return tags + [dv.phenomenon_tag(match.group(1), match.group(2)), dv.TAG_SECTOR_DATA, dv.TAG_PHENOMENA_SET]

    if path == "/sectors":
        return tags + [dv.TAG_SECTOR_DATA, dv.TAG_SYSTEM_SET, dv.TAG_PHENOMENA_SET, dv.TAG_POPULATION, dv.TAG_POSITIONS]
    if path in ("/systems", "/uncharted-systems"):
        return tags + [dv.TAG_SECTOR_DATA, dv.TAG_SYSTEM_DATA, dv.TAG_SYSTEM_ANY, dv.TAG_POPULATION, dv.TAG_POSITIONS]
    if path == "/phenomena":
        return tags + [dv.TAG_SECTOR_DATA, dv.TAG_PHENOMENA_DATA]
    if path.startswith(("/polities", "/species", "/territories", "/planets/")):
        return tags + [dv.TAG_POPULATION, dv.TAG_SYSTEM_DATA]
    return tags + list(ALL_TAGS)


_SECTOR_WRITES = (
    (re.compile(rf"/sectors/{_ID}/(contents|regenerate)"), "sector-contents"),
    (re.compile(rf"/sectors/{_ID}/generate-neighborhood"), "generated"),
    (re.compile(rf"/sectors/{_ID}/wiki"), "sector-edited"),
)


def change_for_write(method, path):
    """
    The change a successful write through the API makes.

    Args:
        method (str): `POST`, `PATCH`, `PUT` or `DELETE`.
        path (str): The path under `/api`.

    Returns:
        tuple: `(kind, ids)`: a key of `dataversion.CHANGES` and the keyword ids `dataversion.record` takes
            (`sector_id`, `system_id`, `phenomenon`), with printed IDs. An unknown write is
            `("everything", {})`.
    """
    for pattern, kind in _SECTOR_WRITES:
        match = pattern.fullmatch(path)
        if match:
            return kind, {"sector_id": match.group(1)}
    if path == "/sectors" and method == "POST":
        return "sector-set", {}
    match = re.fullmatch(rf"/sectors/{_ID}", path)
    if match:
        return ("sector-set" if method == "DELETE" else "sector-edited"), {"sector_id": match.group(1)}

    if path == "/systems" and method == "POST":
        return "system-set", {}
    match = re.fullmatch(rf"/systems/{_ID}(?:/(wiki|star))?", path)
    if match:
        return ("system-set" if method == "DELETE" else "system-edited"), {"system_id": match.group(1)}
    if re.fullmatch(r"/(stars|planets|moons|facilities)(?:/[^/]+(?:/class|/regenerate|/suggest-ratings)?)?", path):
        return "system-any", {}

    match = re.fullmatch(rf"/phenomena/{_ID}/{_ID}(?:/regenerate)?", path)
    if match:
        kind = "phenomena-set" if method == "DELETE" else "phenomenon-edited"
        return kind, {"phenomenon": (match.group(1), match.group(2))}
    return "everything", {}
