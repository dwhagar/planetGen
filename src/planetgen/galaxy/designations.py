"""
planetgen/galaxy/designations.py

The rules of designation history (GEN.1138; docs/design/architecture/designation-history.md): an object may change its name and its ID, the
designations it stopped having are kept, and a retired ID is reserved forever. This module is the pure part: it reads no database and no stream.
The `designations` table and the writers that call these rules follow in the migration of the same item.

An ID is an integer of `object_uid`; a name is text. The retired rows of one designation are `Retired` records.
"""
from typing import NamedTuple

from planetgen.galaxy import object_uid

KIND_UID = "uid"
KIND_NAME = "name"
KINDS = (KIND_UID, KIND_NAME)

REASONS = ("collision bump", "rename", "regeneration", "parent change", "admin renumber")
"""tuple: Why a designation was retired. A moved rogue is not among them: it keeps the ID of its birth position (Boss, 5:49 PM EDT 2026-10-30)."""

CURRENT = "current"
REDIRECT = "redirect"
GONE = "gone"
UNKNOWN = "unknown"
"""str: What a lookup of a designation finds: the live one, a retired one of an object that exists (a redirect), a retired one of a deleted object,
or nothing."""


class DesignationError(ValueError):
    """A change the rules refuse."""


class Retired(NamedTuple):
    """One retired designation. `object_kind` and `object_row_id` name the object that held it (the table and its row number, which never change);
    `gone_epoch` is `None` while the object exists."""

    kind: str
    designation: object
    object_kind: str
    object_row_id: int
    valid_from_epoch: object
    retired_epoch: int
    reason: str
    gone_epoch: object = None


class Resolution(NamedTuple):
    """The answer to a lookup: a `status` and the object it names, with the retired rows (`formerly`) that led to it, newest first."""

    status: str
    object_kind: object = None
    object_row_id: object = None
    formerly: tuple = ()


READ_METHODS = ("GET", "HEAD")
"""tuple: The HTTP methods that only read: they follow a retired ID to the object, where a write refuses it."""

STATUS_REDIRECT = 301
STATUS_CONFLICT = 409
STATUS_GONE = 410


def status_of(method, resolution):
    """
    The HTTP status a request for a designation gets (invariant 6 and the redirect rule): a read of a retired designation of an object that
    exists is a redirect to the object's current one (301); a write against it is refused (409) and does nothing, because a write never acts on a
    guess; the retired designation of a deleted object answers 410 to both. `None` for a current or unknown designation, which the caller's own
    rules (the object, or a 404) answer.
    """
    if resolution.status == REDIRECT:
        return STATUS_REDIRECT if str(method).upper() in READ_METHODS else STATUS_CONFLICT
    if resolution.status == GONE:
        return STATUS_GONE
    return None


def check_reason(reason):
    if reason not in REASONS:
        raise DesignationError(f"a designation is retired for one of {REASONS}, not {reason!r}")
    return reason


def check_new_uid(new_uid, live_uids, retired_uids):
    """
    Refuses an ID that is, or ever was, any object's: invariant 1 of the design. `live_uids` and `retired_uids` are containers of integers (the
    live columns and `designations`' retired IDs); the real check is the unique indexes of the two, and this says why in words first.

    Raises:
        DesignationError: The ID is live or retired.
    """
    if new_uid in live_uids:
        raise DesignationError(f"{object_uid.format_id(new_uid)} is already an object's ID")
    if new_uid in retired_uids:
        raise DesignationError(f"{object_uid.format_id(new_uid)} was an object's ID once and is reserved for ever")
    return new_uid


def cascade(old_system_uid, new_system_uid, body_uids):
    """
    The new IDs of the bodies of a system whose ID changes (invariant 4): a body's ID is its system's with its own body number, so
    `{old body ID: new body ID}` keeps every number. The system itself is not in `body_uids`.

    Raises:
        DesignationError: A body that is not of the old system, the system's own ID in the list, or an unchanged system ID.
    """
    if old_system_uid == new_system_uid:
        raise DesignationError("a renumbering gives the system an ID it does not have")
    if object_uid.body_of(old_system_uid) != 0 or object_uid.body_of(new_system_uid) != 0:
        raise DesignationError("a system's own ID has body number 0")
    moved = {}
    for body in body_uids:
        if object_uid.body_of(body) == 0:
            raise DesignationError(f"{object_uid.format_id(body)} is the system's own ID, not a body's")
        if not object_uid.same_system(body, old_system_uid):
            raise DesignationError(f"{object_uid.format_id(body)} is not a body of {object_uid.format_id(old_system_uid)}")
        moved[body] = object_uid.with_body(new_system_uid, object_uid.body_of(body))
    return moved


def resolve_uid(uid, live_holder, retired_rows):
    """
    What an ID names (invariant 6 and the redirect rule). `live_holder(uid)` gives `(object_kind, row_id)` or `None`; `retired_rows` are the
    `Retired` records of `uid` (at most one: an ID is reserved to one object).

    Returns:
        Resolution: `CURRENT` for a live ID; `REDIRECT` for a retired ID of an object that exists (the caller reads that object's current ID);
            `GONE` for a retired ID of a deleted object; `UNKNOWN` for an ID nobody held.
    """
    holder = live_holder(uid)
    if holder is not None:
        return Resolution(CURRENT, holder[0], holder[1])
    rows = [row for row in retired_rows if row.kind == KIND_UID and row.designation == uid]
    if not rows:
        return Resolution(UNKNOWN)
    row = max(rows, key=lambda r: r.retired_epoch)
    if row.gone_epoch is not None:
        return Resolution(GONE, row.object_kind, row.object_row_id, (row,))
    return Resolution(REDIRECT, row.object_kind, row.object_row_id, (row,))


def resolve_name(name, current_holders, retired_rows):
    """
    What a name names. A name is unique among current holders only (the collision bump hands a stored object's name to another), so the live
    holder wins; failing that the latest past holder that still exists; failing that the latest of a deleted one. `current_holders` is a list of
    `(object_kind, row_id)` (more than one for body names, which are not unique); `retired_rows` the `Retired` records of the name.

    Returns:
        Resolution: with every other past holder in `formerly`, newest first, so a page can offer them.
    """
    if current_holders:
        return Resolution(CURRENT, *current_holders[0])
    rows = sorted((row for row in retired_rows if row.kind == KIND_NAME and row.designation == name),
                  key=lambda r: r.retired_epoch, reverse=True)
    if not rows:
        return Resolution(UNKNOWN)
    alive = [row for row in rows if row.gone_epoch is None]
    chosen = alive[0] if alive else rows[0]
    return Resolution(GONE if not alive else REDIRECT, chosen.object_kind, chosen.object_row_id, tuple(rows))


def resolve_member(text, resolve_parent, current_uid_of):
    """
    What a belt member's or captured body's designation (`<parent ID>.<rank hex>`) names after its parent's ID changed (invariant 5): the parent
    is resolved by `resolve_parent(uid) -> Resolution` and, for a redirect, its current ID read by `current_uid_of(object_kind, row_id)`; the rank is
    kept.

    Returns:
        tuple: `(Resolution of the parent, current designation text or None)`; the text is `None` for an unknown or deleted parent and
            the member's own current text for a live or redirected one.

    Raises:
        ValueError: Text that is not a member designation.
    """
    parent, rank = object_uid.parse_designation(text)
    found = resolve_parent(parent)
    if found.status == CURRENT:
        return found, object_uid.designation(parent, rank)
    if found.status == REDIRECT:
        return found, object_uid.designation(current_uid_of(found.object_kind, found.object_row_id), rank)
    return found, None


def former_lines(rows):
    """
    The designations an object used to have, newest first, as `(kind, designation, reason, retired_epoch)` for the "Formerly" lines of its
    page and the `former_designations` field of its JSON. An ID is printed (`object_uid.format_id`); a name is kept as text.
    """
    out = []
    for row in sorted(rows, key=lambda r: r.retired_epoch, reverse=True):
        text = object_uid.format_id(row.designation) if row.kind == KIND_UID else row.designation
        out.append((row.kind, text, row.reason, row.retired_epoch))
    return out
