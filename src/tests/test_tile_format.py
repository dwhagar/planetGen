# tests/test_tile_format.py

"""
PERF.133: the stored tiles are kept across releases, so a change to what a tile holds or how it is packed has to
say so by raising `planetgen.db.query.TILE_FORMAT` (which rebuilds the stored tiles once). This pins the bytes one fixed
tile packs to for the current format: a change that alters them fails here until `TILE_FORMAT` goes up and the new
digest is recorded beside it.
"""

import hashlib

from planetgen.db import query
from planetgen.web.lib import tilecache

TILE = {
    "placed": [{"id": 1, "x": 1.23456789, "y": 0.0, "z": 0.0}], "planned": [{"shell": 1, "band": 0, "wedge": 2}],
    "filled": {"g": 1}, "generated": [],
    "stars": [{
        "id": 7, "name": "Named", "x": 1.23456789, "y": -2.0, "z": 0.00049, "luminosity_sol": 1234.56789,
        "temperature_k": 5778.4, "radius_sol": 1.234567, "star_type": "G2V Yellow Main Sequence Star",
        "population": "young", "yerkes_class": "V", "shell_index": 3, "band_index": 1, "wedge_index": 4,
        "system_id": None}],
    "points": [], "clouds": [{"id": 3, "x": 1.23456789, "y": 0.0, "z": 0.0, "radius_pc": 2.34567891}],
}

DIGESTS = {
    1: "528fedbb24ba4ff5e32246d8d579b3ef",
}
"""dict: `TILE_FORMAT` to the SHA-256 (first 32 hex) of `TILE`'s stored bytes under it."""


def test_the_stored_bytes_of_a_tile_match_the_recorded_format():
    packed = tilecache._tile_bytes(TILE)
    digest = hashlib.sha256(packed).hexdigest()[:32]
    assert query.TILE_FORMAT in DIGESTS, f"record the digest {digest} for TILE_FORMAT {query.TILE_FORMAT}"
    assert digest == DIGESTS[query.TILE_FORMAT], (
        f"what a tile stores changed (digest {digest}): raise planetgen.db.query.TILE_FORMAT and record "
        f"{digest} for the new number in DIGESTS, so the stored tiles are rebuilt once (PERF.133)")
