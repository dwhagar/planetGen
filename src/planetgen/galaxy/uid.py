# planetgen/galaxy/uid.py

"""
A unique ID for every object in the galaxy (GEN.69, GEN.68's plan).

    sector    its designation (`geometry.provisional_sector_designation`):
              ring, biased layer and slot packed into one integer. It
              needs no row, so a sector nobody has generated has its ID
              already; `sector_uid` is the integer, `format_sector_uid`
              the hex the site prints.
    system, phenomenon
              96 bits (`GALAXY_BITS`): the first 96 bits of
              SHA-256(galaxy seed || "uid-<kind>:" parent "/" index), with
              the top bit set (`derived_uid`) so it can never equal a
              GEN.64 position ID, which fits in 76 bits and which an
              interstellar object or bright-sweep system keeps
              (`names/object_id.py`).
    star, planet, moon, belt, comet
              64 bits (`LOCAL_BITS`) from the same hash, unique under
              their system (a moon under its planet's own ID too), which is
              where they are looked up.

`parent` is the parent's ID in hex and `index` the object's slot in it as
generated (a system's rank in its sector, a planet's orbital index), so
regenerating a unit in place from the same galaxy seed gives every object
the ID it had, and nothing about where an object sits is an input. This is
the same hash the per-unit random seeds use (`galaxy/seed.py`), so there is
nothing new to explain. Nothing here reads a database.
"""

from planetgen.galaxy import seed as galaxy_seed
from planetgen.galaxy.geometry import parse_sector_designation, provisional_sector_designation

GALAXY_BITS = 96
"""int: Width of a system's or phenomenon's ID: unique galaxy-wide."""

LOCAL_BITS = 64
"""int: Width of a star's, planet's, moon's, belt's or comet's ID: unique
under its system."""

NO_SEED = bytes(galaxy_seed.SEED_BYTES)
"""bytes: The seed IDs are drawn from when nothing has been planned (a
one-off system saved on its own)."""


def sector_uid(ring_index, layer_index, slot_index):
    """A sector's ID: its designation as an integer. Raises `ValueError`
    like `provisional_sector_designation` for an address outside the grid."""
    return int(provisional_sector_designation(ring_index, layer_index, slot_index), 16)


def sector_address(uid):
    """The `(ring_index, layer_index, ring_slot_index)` a sector ID names."""
    return parse_sector_designation(f"{uid:X}")


def format_sector_uid(uid):
    """A sector ID as the uppercase hex its designation is printed in."""
    return f"{uid:X}"


def derived_uid(seed, kind, parent, index, bits=GALAXY_BITS):
    """
    One object's ID.

    Args:
        seed (bytes or None): The galaxy's 16-byte seed (`NO_SEED` when
            `None`).
        kind (str): What the object is: "system", "star", "planet",
            "moon", "belt", "comet", or a phenomenon's table.
        parent (str): The parent's ID as text (`format_uid`), or any text
            that names the parent when it has none.
        index (int): The object's slot in its parent.
        bits (int): `GALAXY_BITS` or `LOCAL_BITS`.

    Returns:
        int: `bits` bits; at `GALAXY_BITS` the top bit is set.
    """
    value = galaxy_seed.short_seed(NO_SEED if seed is None else seed, f"uid-{kind}", (parent, index), bits)
    if bits == GALAXY_BITS:
        value |= 1 << (GALAXY_BITS - 1)
    return value


def format_uid(uid, bits=GALAXY_BITS):
    """An ID as fixed-width uppercase hex: 24 digits at `GALAXY_BITS`, 16
    at `LOCAL_BITS`."""
    return f"{uid:0{bits // 4}X}"


def uid_bytes(uid, bits=GALAXY_BITS):
    """An ID as the big-endian bytes a `BINARY(bits/8)` column holds."""
    return uid.to_bytes(bits // 8, "big")


def uid_from_bytes(raw):
    """The integer a `BINARY` column's bytes hold."""
    return int.from_bytes(bytes(raw), "big")
