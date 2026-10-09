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
