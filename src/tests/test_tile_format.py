# tests/test_tile_format.py

"""
PERF.133: the stored nodes are kept across releases, so a change to what a node holds or how it is packed has to
say so by raising `planetgen.db.query.TILE_FORMAT` (which rebuilds the stored nodes once). This pins the bytes one fixed
node packs to for the current format: a change that alters them fails here until `TILE_FORMAT` goes up and the new
digest is recorded beside it.
"""

import hashlib

from planetgen.db import query
from planetgen.web.lib import tilecache

NODE = {
    "placed": [{"id": 1, "x": 1.23456789, "y": 0.0, "z": 0.0}], "planned": [{"shell": 1, "band": 0, "wedge": 2}],
    "filled": {"g": 1}, "generated": [],
    "aggregate": {"generated": 5, "stats": {"systems": 3.0}, "lights": {"stars": 4, "star_luminosity_sol": 1234.56789}},
    "stars": [{
        "id": 7, "name": "Named", "x": 1.23456789, "y": -2.0, "z": 0.00049, "luminosity_sol": 1234.56789,
        "temperature_k": 5778.4, "radius_sol": 1.234567, "star_type": "G2V Yellow Main Sequence Star",
        "population": "young", "yerkes_class": "V", "shell_index": 3, "band_index": 1, "wedge_index": 4,
        "system_id": None}],
    "points": [], "clouds": [{"id": 3, "x": 1.23456789, "y": 0.0, "z": 0.0, "radius_pc": 2.34567891}],
}

DIGESTS = {
    3: "365804129a9ae78d5135a26d121616ed",
    4: "357d5259c09ea2060f7b9d665b6e64ef",
}
"""dict: `TILE_FORMAT` to the SHA-256 (first 32 hex) of `NODE`'s stored bytes under it."""


def test_the_stored_bytes_of_a_node_match_the_recorded_format():
    packed = tilecache._node_bytes(NODE)
    digest = hashlib.sha256(packed).hexdigest()[:32]
    assert query.TILE_FORMAT in DIGESTS, f"record the digest {digest} for TILE_FORMAT {query.TILE_FORMAT}"
    assert digest == DIGESTS[query.TILE_FORMAT], (
        f"what a node stores changed (digest {digest}): raise planetgen.db.query.TILE_FORMAT and record "
        f"{digest} for the new number in DIGESTS, so the stored nodes are rebuilt once (PERF.133)")
