# planetgen/names/naming_key.py

"""
The naming key (GEN.70; docs/design/object-ids.md, "Planned" step 3).

A galaxy has one naming key: 8 uppercase hex digits kept in the control
database (`galaxy_naming`, one row per galaxy database). It is drawn from
the galaxy seed when the galaxy is planned and an admin can change it.
The objects with no star-derived name (rogue planets, standalone black
holes and neutron stars, nebulae, supernova remnants and their cores,
quasars, interstellar comets, asteroid fields, constellations) are named
by `codec_name`: the phoneme codec applied to the object's ID in a domain
made of the key and the object's kind. Nothing about a name is stored, so
changing the key renames all of them at once. Stars, sectors, planets,
moons and belts do not use the key.

The codec turns its domain into an affine permutation of the hex digits
(`GatedPhonemeCodec._get_domain_params`), so a key moves a name among 128
permutations per kind; two different keys can therefore name an object
alike. The ID stays the identity.

`CODEC_VERSION` goes up whenever the codec would spell an ID differently
(its golden words in `test_gated_phoneme_codec.py` change); it is stored
beside the key so an admin sees when names were made under another one.
"""

import datetime
import hashlib
import re

from planetgen.names import object_id
from planetgen.names.gated_phoneme_codec import GatedPhonemeCodec

CODEC_VERSION = 1
"""int: The codec's version; see the module docstring."""

KEY_DIGITS = 8
"""int: A naming key's length in hex digits (32 bits)."""

_KEY = re.compile(r"[0-9A-Fa-f]{%d}" % KEY_DIGITS)

OBJECT_ID_DIGITS = 19
"""int: Every GEN.64 object ID is 19 hex digits; the decoder is told so
(`decode(phrase, domain, length=19)` is exact)."""

KINDS = (
    "rogue-planet", "black-hole", "neutron-star", "nebula", "supernova-remnant",
    "quasar", "comet", "asteroid-field", "black-hole-core", "neutron-star-core", "constellation",
)
"""tuple: The kinds the codec names: `object_id.KIND_CODES`' spelling
(less the bright-sweep system, which is a star system), and
`"constellation"` (VIEW.4)."""

_CODEC = GatedPhonemeCodec()


def draw_key(galaxy_seed):
    """A galaxy's first naming key: the first 32 bits of
    SHA-256(galaxy seed || "naming-key"), so the same seed always gets the
    same key. `galaxy_seed` is the 16-byte seed."""
    return hashlib.sha256(bytes(galaxy_seed) + b"naming-key").hexdigest()[:KEY_DIGITS].upper()


def parse_key(text):
    """
    A key typed as 8 hex digits (either case), returned in upper case.

    Raises:
        ValueError: If `text` isn't exactly 8 hex digits.
    """
    text = str(text).strip()
    if not _KEY.fullmatch(text):
        raise ValueError(f"a naming key is {KEY_DIGITS} hex digits (0-9, A-F), not {text!r}")
    return text.upper()


def codec_domain(key, kind):
    """The codec domain for objects of `kind` under `key`."""
    if kind not in KINDS:
        raise ValueError(f"no codec name for the kind {kind!r}")
    return f"{parse_key(key)}:{kind}"


def codec_name(object_id, kind, key):
    """
    The name of the object with 19-digit hex `object_id` (`object_id.format_id`)
    of `kind`, under naming key `key`: capitalized words, e.g. "Bafor Telun".
    """
    if len(object_id) != OBJECT_ID_DIGITS:
        raise ValueError(f"an object ID is {OBJECT_ID_DIGITS} hex digits, not {object_id!r}")
    return _CODEC.encode(object_id, domain=codec_domain(key, kind)).title()


def object_id_of(name, kind, key):
    """The 19-digit hex ID a `codec_name` stands for (its inverse)."""
    return _CODEC.decode(name.lower(), codec_domain(key, kind), length=OBJECT_ID_DIGITS)


# ---------------------------------------------------------------------
# Names shown for stored object IDs (GEN.71)
# ---------------------------------------------------------------------

_resolver = None


def set_resolver(function):
    """Installs the callable that answers "what is the naming key in force
    for this request?" (the API sets one that reads the control database;
    without one, or when it answers `None`, an ID is shown as the ID)."""
    global _resolver
    _resolver = function


def active_key():
    """The naming key in force, or `None` (no resolver, or none drawn)."""
    return _resolver() if _resolver is not None else None


def codec_kind_of(stored_name):
    """
    The kind (`KINDS`) the 19-hex-digit object ID `stored_name` belongs to,
    or `None` when it isn't an object ID or is one of a kind the codec
    doesn't name (a bright-sweep system keeps its ID as its uid).
    """
    if not isinstance(stored_name, str) or len(stored_name) != OBJECT_ID_DIGITS:
        return None
    try:
        kind = object_id.parse_id(stored_name)["kind"]
    except ValueError:
        return None
    return kind if kind in KINDS else None


def display_name(stored_name, key=None):
    """
    What to show for `stored_name`: the codec name of an object ID under
    `key` (default: `active_key()`), anything else unchanged. The ID stays
    stored; the kind comes from the ID's own type bits.
    """
    kind = codec_kind_of(stored_name)
    if kind is None:
        return stored_name
    key = key if key is not None else active_key()
    if key is None:
        return stored_name
    return codec_name(stored_name.upper(), kind, key)


def rename_names(value, get_key):
    """
    `value` (JSON-able dicts and lists) with every `"name"` that is an
    object ID replaced by its codec name. `get_key` is called at most once,
    and only when such a name turns up, so payloads without one cost no
    key lookup. Returns `value` itself when nothing was renamed.
    """
    state = {}

    def key():
        if "key" not in state:
            state["key"] = get_key()
        return state["key"]

    def walk(item):
        if isinstance(item, dict):
            changed = None
            for field, child in item.items():
                if field == "name" and codec_kind_of(child) is not None:
                    renamed = display_name(child, key())
                else:
                    renamed = walk(child) if isinstance(child, (dict, list)) else child
                if renamed is not child:
                    if changed is None:
                        changed = dict(item)
                    changed[field] = renamed
            return item if changed is None else changed
        if isinstance(item, list):
            out = None
            for index, child in enumerate(item):
                renamed = walk(child) if isinstance(child, (dict, list)) else child
                if renamed is not child:
                    if out is None:
                        out = list(item)
                    out[index] = renamed
            return item if out is None else out
        return item

    return walk(value)


def stored_name_for(typed, key=None):
    """
    The object ID a typed codec name stands for under `key` (default:
    `active_key()`), or `None`: the name is tried in each codec kind's
    domain and kept when the decoded ID is of that kind.
    """
    key = key if key is not None else active_key()
    if key is None or not isinstance(typed, str) or not typed.strip():
        return None
    for kind in KINDS:
        try:
            decoded = object_id_of(" ".join(typed.lower().split()), kind, key)
        except Exception:  # noqa: BLE001 -- the codec rejects what it can't read
            continue
        if codec_kind_of(decoded) == kind:
            return decoded
    return None


# ---------------------------------------------------------------------
# The control database's row (one per galaxy database)
# ---------------------------------------------------------------------

def get(conn, database):
    """`{"key", "codec_version", "drawn_at", "changed_at", "changed_by"}` for
    `database` from an open control connection, or `None` before one is drawn."""
    row = conn.execute(
        "SELECT naming_key, codec_version, drawn_at, changed_at, changed_by FROM galaxy_naming"
        " WHERE database_name = ?", (database,)).fetchone()
    if row is None:
        return None
    return {"key": row["naming_key"], "codec_version": row["codec_version"], "drawn_at": row["drawn_at"],
            "changed_at": row["changed_at"], "changed_by": row["changed_by"]}


def key_of(conn, database):
    """The key in force for `database`, or `None` before one is drawn."""
    stored = get(conn, database)
    return stored["key"] if stored else None


def draw(conn, database, galaxy_seed, replace=False):
    """
    Stores the key drawn from `galaxy_seed` for `database` when it has none
    (or when `replace`, for a new galaxy planned over an old one), and
    returns the key in force.
    """
    existing = get(conn, database)
    if existing is not None and not replace:
        return existing["key"]
    key = draw_key(galaxy_seed)
    conn.execute(
        "INSERT INTO galaxy_naming (database_name, naming_key, codec_version, drawn_at)"
        " VALUES (?, ?, ?, NOW(6)) ON DUPLICATE KEY UPDATE naming_key = ?,"
        " codec_version = ?, drawn_at = NOW(6), changed_at = NULL, changed_by = NULL",
        (database, key, CODEC_VERSION, key, CODEC_VERSION))
    conn.commit()
    return key


def change(conn, database, key, changed_by):
    """
    Sets `database`'s naming key to `key` (8 hex digits) on behalf of the
    admin `changed_by`, and returns the key stored.

    Raises:
        ValueError: If `key` isn't 8 hex digits.
        LookupError: If no key was drawn yet (the galaxy isn't planned).
    """
    key = parse_key(key)
    if get(conn, database) is None:
        raise LookupError(f"{database} has no naming key yet; plan the galaxy first")
    conn.execute(
        "UPDATE galaxy_naming SET naming_key = ?, codec_version = ?, changed_at = ?, changed_by = ?"
        " WHERE database_name = ?",
        (key, CODEC_VERSION, datetime.datetime.now(datetime.timezone.utc).replace(tzinfo=None),
         str(changed_by)[:64], database))
    conn.commit()
    return key
