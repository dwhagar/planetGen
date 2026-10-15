# tests/maps_pending.py

"""The marker for tests of the 3D Galaxy Map and its tile, box and stage reads: they wait for Maps to rebuild
them on shells, bands and wedges (MAP.205/206), and Maps removes the marker as it restores each one (UX.111)."""

import pytest

MAPS_PENDING = pytest.mark.skip(
    reason="UX.111: the 3D Galaxy Map is rebuilt for shells, bands and wedges by Maps (MAP.205/206); "
           "Maps restores these with it")
