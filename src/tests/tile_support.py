"""Reading a `/galaxy/tiles` answer in tests (MAP.159: it is packed, not JSON)."""

from planetgen.web.lib import tilepack


def tile_payload(response):
    """What `tilecache.fetch_tiles` returns, from a test client's `/galaxy/tiles` response:
    the header's fields and `tiles` as dicts."""
    assert response.mimetype == tilepack.MIMETYPE, response.get_data(as_text=True)[:200]
    header, tiles = tilepack.unpack_response(response.get_data())
    return {**header, "tiles": tiles}
