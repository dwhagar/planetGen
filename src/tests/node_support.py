"""Reading a `/galaxy/nodes` answer in tests (MAP.159: it is packed, not JSON)."""

from planetgen.web.lib import tilepack


def node_payload(response):
    """What `tilecache.fetch_nodes` returns, from a test client's `/galaxy/nodes` response:
    the header's fields and `nodes` as dicts."""
    assert response.mimetype == tilepack.MIMETYPE, response.get_data(as_text=True)[:200]
    header, nodes = tilepack.unpack_response(response.get_data())
    return {**header, "nodes": nodes}
