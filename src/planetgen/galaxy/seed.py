# planetgen/galaxy/seed.py

"""
The galaxy's 128-bit seed and the seeds derived from it (GEN.39; see
docs/design/generation/reproducible-galaxies.md, section 3).

One seed belongs to the galaxy: 16 bytes in `galaxy_shape.galaxy_seed`
(`BINARY(16)`), shown and typed as 32 hex digits, written once when the
galaxy is first planned (`planetgen plan --seed`, or drawn at random).
Every unit of work gets its own seed from it:

    unit seed = SHA-256(galaxy seed || "kind:" address)

for example `sector:12/3/0` (a sector's address is shell/band/wedge). The whole 256-bit digest seeds the unit's
own draw stream (`seeded`, see `planetgen.util.draw`), so a sector's numbers depend
only on the galaxy seed and its address, never on the order sectors run
in or how many workers run them. The version isn't mixed in: a release
that changes one formula changes only what that formula touches.
"""

import contextlib
import contextvars
import hashlib
import re
import secrets

from planetgen.util import draw


SEED_BYTES = 16
"""int: A galaxy seed's length: 128 bits (Boss, 2026-10-02)."""

SEED_HEX_DIGITS = 2 * SEED_BYTES
"""int: A galaxy seed as text: 32 hex digits."""

_HEX = re.compile(r"[0-9A-Fa-f]{%d}" % SEED_HEX_DIGITS)

def new_seed():
    """A fresh galaxy seed: 16 bytes from the operating system's random
    source, the only draw that may come from there."""
    return secrets.token_bytes(SEED_BYTES)


def parse_seed(text):
    """
    A galaxy seed typed as 32 hex digits (either case; spaces and a `0x`
    prefix are not accepted).

    Returns:
        bytes: The 16-byte seed.

    Raises:
        ValueError: If `text` isn't exactly 32 hex digits.
    """
    text = str(text).strip()
    if not _HEX.fullmatch(text):
        raise ValueError(f"a galaxy seed is {SEED_HEX_DIGITS} hex digits (0-9, A-F), not {text!r}")
    return bytes.fromhex(text)


def format_seed(seed):
    """A 16-byte galaxy seed as 32 uppercase hex digits."""
    return bytes(seed).hex().upper()


def address_text(address):
    """A unit's address as its seed text spells it: a tuple's parts
    joined by `/` (`(12, 3, 0)` -> `12/3/0`, shell 12, band 3, wedge 0), anything else as `str`."""
    if isinstance(address, (tuple, list)):
        return "/".join(str(part) for part in address)
    return str(address)


def unit_seed(galaxy_seed, kind, address):
    """
    One unit's seed: the SHA-256 of the galaxy seed's 16 bytes followed by
    `"kind:address"` in UTF-8, as a 256-bit integer.

    Args:
        galaxy_seed (bytes): The galaxy's 16-byte seed.
        kind (str): The unit's kind, for example `"sector"`.
        address: Its address (see `address_text`).

    Returns:
        int: The digest, big-endian.
    """
    if len(galaxy_seed) != SEED_BYTES:
        raise ValueError(f"a galaxy seed is {SEED_BYTES} bytes, not {len(galaxy_seed)}")
    digest = hashlib.sha256(bytes(galaxy_seed) + f"{kind}:{address_text(address)}".encode("utf-8")).digest()
    return int.from_bytes(digest, "big")


def short_seed(galaxy_seed, kind, address, bits=63):
    """The top `bits` bits of `unit_seed`, for a seed stored in a narrower
    column (the bright-star scatter's `BIGINT UNSIGNED`)."""
    return unit_seed(galaxy_seed, kind, address) >> (256 - bits)


_current_galaxy_seed = contextvars.ContextVar("galaxy_seed", default=None)


def current_galaxy_seed():
    """The galaxy seed of the unit running in this thread (`seeded` binds it), or `None` outside any unit: work with
    no galaxy, a stand-alone command or a test. Hill spheres read the galaxy's mass chain through it (GEN.400)."""
    return _current_galaxy_seed.get()


_current_disk_scale_length_pc = contextvars.ContextVar("disk_scale_length_pc", default=None)


def current_disk_scale_length_pc():
    """The thin disk's scale length, parsecs, of the galaxy the running unit belongs to, when `seeded` was given it
    (a planned galaxy's real shape, whatever options drew it), or `None`: `galaxy.scales` then takes the length the seed
    draws."""
    return _current_disk_scale_length_pc.get()


_current_galaxy_shape = contextvars.ContextVar("galaxy_shape", default=None)


def current_galaxy_shape():
    """The stored `GalaxyShape` of the galaxy the running unit belongs to, when `seeded` was given it, or `None`: the
    Hill spheres' enclosed mass (`galaxy.enclosed_mass`) then uses the shape the galaxy seed draws."""
    return _current_galaxy_shape.get()


@contextlib.contextmanager
def seeded(galaxy_seed, kind, address, disk_scale_length_pc=None, shape=None):
    """
    Runs one unit on its own seed: binds a `draw.Stream` seeded from
    `unit_seed` for the length of the block (in this thread only), so
    what ran before and after carries on as if the unit hadn't run here.
    With no galaxy seed (`None`: nothing planned) it changes nothing.

    Args:
        disk_scale_length_pc (float, optional): The galaxy's thin disk scale length, which the laws quoted in the
            galaxy's own scales (GEN.404) read through `galaxy.scales`.
        shape (GalaxyShape, optional): The galaxy's stored shape, which the Hill spheres' enclosed mass reads through
            `galaxy.enclosed_mass` (GEN.370); without it the shape the seed draws stands in.
    """
    if galaxy_seed is None:
        yield
        return
    token = _current_galaxy_seed.set(bytes(galaxy_seed))
    length_token = _current_disk_scale_length_pc.set(disk_scale_length_pc)
    shape_token = _current_galaxy_shape.set(shape)
    try:
        with draw.bound(unit_seed(galaxy_seed, kind, address)):
            yield
    finally:
        _current_galaxy_shape.reset(shape_token)
        _current_disk_scale_length_pc.reset(length_token)
        _current_galaxy_seed.reset(token)
