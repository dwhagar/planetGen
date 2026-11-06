# planetgen/galaxy/gravity_cache.py

"""
The gravity map's grid cache: file keys, stamps and the read and write of one sector's grid (MAP.387; docs/design/orbits/gravity-map.md)
=========================================================================================================================================

A sector's grid (`db.gravity_grid`) is derived data: it is worked out from the stored stars and thrown away when they move, so it is never a galaxy table and
needs no migration. It is kept as one file per sector in a cache directory the caller names. This module is the part the web layer may read and does not need
numba: the key of a file, the stamp that makes the key change, and the file itself.

- **Stamp.** `stamp(parts)` is 16 hex characters of the SHA-256 of `parts`, which the job builds from the rule version, the grid's settings, the potential and, for
  every sector in the grid's reach, the sector's `modified_at` and its rows' count and latest `modified_at`. The orbital update (GEN.105) writes the rows it moves,
  so a moved star changes the stamp of every grid whose reach holds its sector and of no other; a sector farther than the reach changes nothing here.
- **File.** `g<sector id>_<zones>_<stamp>.bin`: the magic `PGG1`, the zone count and the column count as two unsigned 32-bit integers, then the values as float32,
  `zones^3` rows of `COLUMNS` in `gravity_kernel`'s column order (the pull, the potential and the tidal tensor of the *local* field: the generated stars less the
  smooth galaxy's own stars, MAP.383), zone `(i, j, k)` at row `(i * zones + j) * zones + k`. The reader adds the smooth galaxy to get the total.
- **Failing open.** A missing, short, torn or foreign file reads as `None` (the grid is worked out again); an unwritable directory loses the cache and nothing
  else. A write goes to a temporary file and is renamed, so a reader never sees half a file. Writing a sector's grid removes that sector's older files.
"""

import hashlib
import os
import re
import struct
import tempfile

import numpy as np

MAGIC = b"PGG1"
"""bytes: The first four bytes of a grid file."""
RULE_VERSION = 1
"""int: The version of the rule the grid is worked out by (MAP.383); a change to it changes every stamp."""
COLUMNS = 10
"""int: The columns of a grid row (`galaxy.gravity_kernel`'s)."""
_HEADER = struct.Struct("<4sII")
_NAME = re.compile(r"^g(\d+)_(\d+)_([0-9a-f]{16})\.bin$")


def stamp(parts):
    """16 hex characters naming `parts` (anything with a stable `repr`: numbers, strings, tuples of them)."""
    return hashlib.sha256(repr(parts).encode("utf-8")).hexdigest()[:16]


def file_name(sector_id, zones, stamp_hex):
    """The file name of a sector's grid."""
    return f"g{int(sector_id)}_{int(zones)}_{stamp_hex}.bin"


def read(directory, sector_id, zones, stamp_hex):
    """The grid, a float32 array `(zones^3, COLUMNS)`, or `None` when the file is missing or not a whole grid of this shape."""
    try:
        with open(os.path.join(directory, file_name(sector_id, zones, stamp_hex)), "rb") as handle:
            raw = handle.read()
        magic, found_zones, columns = _HEADER.unpack_from(raw)
        count = found_zones ** 3 * columns
        if magic != MAGIC or found_zones != zones or columns != COLUMNS or len(raw) != _HEADER.size + 4 * count:
            return None
        return np.frombuffer(raw, dtype="<f4", count=count, offset=_HEADER.size).reshape(found_zones ** 3, columns).astype(np.float32)
    except (OSError, struct.error, ValueError):
        return None


def write(directory, sector_id, zones, stamp_hex, values):
    """
    Stores a grid (`(zones^3, COLUMNS)`) and removes the sector's older files. Returns whether it was stored; an unwritable directory is `False`, not an error.
    """
    values = np.ascontiguousarray(values, dtype="<f4")
    if values.shape != (zones ** 3, COLUMNS):
        raise ValueError(f"a grid of {zones} zones is ({zones ** 3}, {COLUMNS}), not {values.shape}")
    try:
        os.makedirs(directory, exist_ok=True)
        handle, temporary = tempfile.mkstemp(dir=directory, prefix=".g", suffix=".tmp")
        try:
            with os.fdopen(handle, "wb") as out:
                out.write(_HEADER.pack(MAGIC, zones, COLUMNS))
                out.write(values.tobytes())
            os.replace(temporary, os.path.join(directory, file_name(sector_id, zones, stamp_hex)))
        except BaseException:
            try:
                os.unlink(temporary)
            except OSError:
                pass
            raise
    except OSError:
        return False
    prune(directory, sector_id, keep=file_name(sector_id, zones, stamp_hex))
    return True


def prune(directory, sector_id, keep=None):
    """Removes every grid file of the sector except `keep`; returns how many went."""
    removed = 0
    try:
        names = os.listdir(directory)
    except OSError:
        return 0
    for name in names:
        found = _NAME.match(name)
        if found and int(found.group(1)) == int(sector_id) and name != keep:
            try:
                os.unlink(os.path.join(directory, name))
                removed += 1
            except OSError:
                pass
    return removed


def clear(directory):
    """Removes every grid file (and stray temporaries) in `directory`; returns how many went."""
    removed = 0
    try:
        names = os.listdir(directory)
    except OSError:
        return 0
    for name in names:
        if _NAME.match(name) or (name.startswith(".g") and name.endswith(".tmp")):
            try:
                os.unlink(os.path.join(directory, name))
                removed += 1
            except OSError:
                pass
    return removed
