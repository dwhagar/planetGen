# planetgen/names/naming_key.py

"""
The naming key (GEN.70; docs/design/architecture/object-ids.md, "Planned" step 3).

A galaxy has one naming key: 8 uppercase hex digits kept in the control
database (`galaxy_naming`, one row per galaxy database). It is drawn from
the galaxy seed when the galaxy is planned and an admin can change it.
The objects with no star-derived name (rogue planets, standalone black
holes and neutron stars, nebulae, supernova remnants and their cores,
quasars, interstellar comets and asteroid fields, GEN.1256) keep their
Object ID as the stored `name` (the printed form, `object_uid.format_id`)
and are shown by `local_name`: two capitalized words spelling the serial and
body number of the ID, scrambled by a keyed bijection of the birth sector.
Nothing about a name is stored, so changing the key renames all of them at
once. Regions of the galaxy are named from their block (`galaxy.named_regions`)
under the same key. Stars, sectors, planets, moons and belts do not use the key.

The codec turns its domain into an affine permutation of the hex digits
(`GatedPhonemeCodec._get_domain_params`), so a key moves a name among 128
permutations per domain; two different keys can therefore name an object
alike. The ID stays the identity.

`CODEC_VERSION` goes up whenever the codec would spell an ID differently
(its golden words in `test_gated_phoneme_codec.py` change); it is stored
beside the key so an admin sees when names were made under another one.
"""

import datetime
import functools
import hashlib
import re

from planetgen.galaxy import object_uid
from planetgen.names.gated_phoneme_codec import GatedPhonemeCodec

CODEC_VERSION = 2
"""int: The codec's version; see the module docstring. 2: local names from the Object ID (GEN.1255)."""

KEY_DIGITS = 8
"""int: A naming key's length in hex digits (32 bits)."""

_KEY = re.compile(r"[0-9A-Fa-f]{%d}" % KEY_DIGITS)

KINDS = ("constellation", "region")
"""tuple: The kinds `codec_domain` names: `"constellation"` (VIEW.4) and `"region"` (GEN.825: a region of the galaxy is named from
its block, `galaxy.named_regions`)."""

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


# ---------------------------------------------------------------------
# Local names from the Object ID (GEN.1255)
# ---------------------------------------------------------------------
#
# A local name spells the serial and body number of an Object ID (`galaxy.object_uid`: 28 + 12 bits, 10 hex digits) and nothing of
# the birth sector, as two words of about 12 letters ("Banhunm Finkap"). Inside a sector it is exact (the spelling is a bijection of the
# 40 bits); on its own it can be shared with an object of another sector, so the sector is shown beside it
# (`local_name_with_sector`) and the object ID stays the identity.
#
# Spelled directly, serials 0, 1, 2 ... share their first word, so the 26 count bits of the serial are first scrambled by a keyed
# bijection, `y = (a * x + b) mod 2**26` with `a` odd, both drawn from the naming key and the sector. `a` is also required to
# keep any two of the first `SPREAD_SERIALS` serials at least `SPREAD_GAP` apart after the scramble, which is what makes their
# first words differ (a first word holds all but the lowest eight bits of the scramble).

LOCAL_DOMAIN = "local"
"""str: The codec domain's kind of a local name; it is not one of `KINDS` (no object ID of that kind exists)."""

LOCAL_DIGITS = (object_uid.LAYOUT.serial_bits + object_uid.LAYOUT.body_bits) // 4
"""int: Hex digits a local name spells: the serial (7) and the body number (3)."""

SPREAD_SERIALS = 4096
"""int: How many of a sector's first serials are guaranteed first words of their own."""

SPREAD_GAP = 256
"""int: The least distance, after the scramble, between any two of the first `SPREAD_SERIALS` serials (the first word drops the low 8 bits)."""

_SCRAMBLE_BITS = object_uid.LAYOUT.serial_count_bits
_SCRAMBLE_MASK = (1 << _SCRAMBLE_BITS) - 1


def local_domain(key):
    """The codec domain of local names under `key`."""
    return f"{parse_key(key)}:{LOCAL_DOMAIN}"


@functools.lru_cache(maxsize=4096)
def _scramble_params(key, sector):
    """The multiplier `a` (odd) and offset `b` of the scramble for `sector` (its 56-bit address as an integer) under `key`: the first
    multiplier drawn from SHA-256(key, sector, attempt) that keeps the first `SPREAD_SERIALS` serials `SPREAD_GAP` apart."""
    sector_text = f"{sector:0{object_uid.LAYOUT.sector_bits // 4}X}"
    attempt = 0
    while True:
        digest = hashlib.sha256(f"local-name:{parse_key(key)}:{sector_text}:{attempt}".encode()).digest()
        a = int.from_bytes(digest[:8], "big") & _SCRAMBLE_MASK | 1
        b = int.from_bytes(digest[8:16], "big") & _SCRAMBLE_MASK
        spread = min(min(value, (1 << _SCRAMBLE_BITS) - value)
                     for value in ((a * step) & _SCRAMBLE_MASK for step in range(1, SPREAD_SERIALS + 1)))
        if spread >= SPREAD_GAP:
            return a, b
        attempt += 1


def scramble(count, sector, key):
    """The count of a serial after the keyed scramble of its sector (`sector` is the 56-bit address as an integer)."""
    a, b = _scramble_params(parse_key(key), sector)
    return (a * count + b) & _SCRAMBLE_MASK


def unscramble(scrambled, sector, key):
    """The count a `scramble` result came from (its inverse)."""
    a, b = _scramble_params(parse_key(key), sector)
    return (pow(a, -1, 1 << _SCRAMBLE_BITS) * (scrambled - b)) & _SCRAMBLE_MASK


def _sector_of(uid):
    layout = object_uid.LAYOUT
    return uid >> (layout.serial_bits + layout.body_bits)


def sector_text(uid):
    """The birth sector of an Object ID (an integer) as it is shown beside a local name: its 14 hex digits, `00FA0FFFE9E240`."""
    return f"{_sector_of(uid):0{object_uid.LAYOUT.sector_bits // 4}X}"


def local_name(uid, key):
    """
    The local name of the object whose Object ID (`galaxy.object_uid`, an integer) is `uid`, under naming key `key`: two capitalized
    words spelling its serial and body number.
    """
    layout = object_uid.LAYOUT
    fields = object_uid.unpack(uid, layout)
    count = scramble(fields.serial, _sector_of(uid), key)
    number = (fields.serial_kind << layout.serial_count_bits) | count
    spelled = f"{(number << layout.body_bits) | fields.body:0{LOCAL_DIGITS}X}"
    return _CODEC.encode(spelled, domain=local_domain(key)).title()


def local_name_with_sector(uid, key):
    """The local name followed by its sector, the way a name that stands alone is shown: `Banhunm Finkap, sector 00FA0FFFE9E240`."""
    return f"{local_name(uid, key)}, sector {sector_text(uid)}"


def uid_of_local_name(name, sector, key):
    """
    The Object ID (an integer) a local name stands for in `sector`, under `key` (the inverse of `local_name`). `sector` is the birth
    sector's 56-bit address as an integer or as its 14 hex digits.

    Raises:
        ValueError: A name the codec cannot read as a local name, or a sector that is not 56 bits.
    """
    layout = object_uid.LAYOUT
    if isinstance(sector, str):
        digits = sector.strip().upper().removeprefix("SECTOR").strip()
        if not re.fullmatch(r"[0-9A-F]{%d}" % (layout.sector_bits // 4), digits):
            raise ValueError(f"a sector is {layout.sector_bits // 4} hex digits, not {sector!r}")
        sector = int(digits, 16)
    if not isinstance(sector, int) or not 0 <= sector < (1 << layout.sector_bits):
        raise ValueError(f"a sector address is {layout.sector_bits} bits, not {sector!r}")
    spelled = int(_CODEC.decode(" ".join(name.lower().split()), local_domain(key), length=LOCAL_DIGITS), 16)
    body = spelled & ((1 << layout.body_bits) - 1)
    number = spelled >> layout.body_bits
    kind = number >> layout.serial_count_bits
    count = unscramble(number & _SCRAMBLE_MASK, sector, key)
    serial = (kind << layout.serial_count_bits) | count
    return (((sector << layout.serial_bits) | serial) << layout.body_bits) | body


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


def stored_uid_of(stored_name):
    """
    The Object ID (an integer) a stored `name` is, or `None`: the objects named by their ID (GEN.1256) keep its printed form
    (`object_uid.format_id`, `FE81000A2B-0000005-000`) in `name`, and no other name has that shape.
    """
    if not isinstance(stored_name, str) or not _STORED_ID.fullmatch(stored_name):
        return None
    try:
        return object_uid.parse_id(stored_name)
    except ValueError:
        return None


_STORED_ID = re.compile(r"[0-9A-Fa-f]+(-[0-9A-Fa-f]+){2}")


def display_name(stored_name, key=None):
    """
    What to show for `stored_name`: the local name of an Object ID under `key` (default: `active_key()`), anything else unchanged.
    The ID stays stored.
    """
    uid = stored_uid_of(stored_name)
    key = key if key is not None else active_key()
    if uid is None or key is None:
        return stored_name
    return local_name(uid, key)


SECTOR_FIELD = "name_sector"
"""str: The field `rename_names` adds beside a renamed `name`: the birth sector's 14 hex digits (`sector_text`)."""


def rename_names(value, get_key):
    """
    `value` (JSON-able dicts and lists) with every `"name"` that is a stored Object ID replaced by its local name, and `SECTOR_FIELD`
    added beside it. `get_key` is called at most once, and only when such a name turns up, so payloads without one cost no key
    lookup. Returns `value` itself when nothing was renamed.
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
                uid = stored_uid_of(child) if field == "name" else None
                if uid is not None and key() is not None:
                    if changed is None:
                        changed = dict(item)
                    changed[field] = local_name(uid, key())
                    changed[SECTOR_FIELD] = sector_text(uid)
                    continue
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


def uid_of_typed_name(typed, key=None):
    """
    The Object ID (an integer) a typed name stands for, or `None`: the name is `Banhunm Finkap, sector 00FA0FFFE9E240`
    (`local_name_with_sector`), the sector after a comma, with or without the word `sector`. Needs the key in force.
    """
    key = key if key is not None else active_key()
    if key is None or not isinstance(typed, str) or "," not in typed:
        return None
    name, _comma, sector = typed.partition(",")
    try:
        return uid_of_local_name(name, sector, key)
    except Exception:  # noqa: BLE001 -- the codec rejects what it can't read
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
