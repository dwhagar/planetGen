# wikiClient/client.py

"""
`WikiClient`: the single object callers construct once per wiki target and
call `create_page` on -- a thin dispatcher over whichever backend
`backend=` names (`"wikijs"` -> `wikijs.WikiJsBackend`, `"mediawiki"` ->
`mediawiki.MediaWikiBackend`). Both backends implement the same
`base.WikiBackend.create_page` contract and raise the same `exceptions`
hierarchy, so calling code never needs to know or care which wiki software
is actually on the other end of a given `WikiClient` -- only which one it
picked when constructing it.

Kept single-backend per instance rather than "publish to both at once":
a system destined for both a Wiki.js and a MediaWiki instance (this
project's schema tracks a separate `wikijs_url`/`mediawiki_url` per
system) just gets two `WikiClient`s, one per target, each fed the matching
pre-rendered content (`markdown_content` for `"wikijs"`, `wikitext_content`
for `"mediawiki"` -- see `wikijs.py`/`mediawiki.py`'s own `create_page`
docstrings).

    from wikiClient import WikiClient

    conduit = WikiClient(backend="wikijs", base_url="https://wiki.example.com", api_token="...")
    page = conduit.create_page(path="systems/kepler-442", title="Kepler-442", content=markdown_content)

    conduit = WikiClient(
        backend="mediawiki", base_url="https://wiki.example.com/w",
        username="Bot@planetgen", password="...",
    )
    page = conduit.create_page(path="Systems/Kepler-442", title="Kepler-442", content=wikitext_content)
"""

from .mediawiki import MediaWikiBackend
from .wikijs import WikiJsBackend

_BACKENDS = {
    "wikijs": WikiJsBackend,
    "mediawiki": MediaWikiBackend,
}


class WikiClient:
    """Constructed with `backend` plus that backend's own connection
    keyword arguments (see `WikiJsBackend`/`MediaWikiBackend`'s own
    `__init__` docstrings for what each expects); every other method call
    is forwarded to that backend unchanged.

    Args:
        backend (str): `"wikijs"` or `"mediawiki"`.
        **kwargs: Forwarded to the chosen backend's constructor.

    Raises:
        ValueError: `backend` isn't one of the known backend names.
    """

    def __init__(self, backend, **kwargs):
        try:
            backend_cls = _BACKENDS[backend]
        except KeyError:
            raise ValueError(f"Unknown wiki backend {backend!r} -- choose one of {sorted(_BACKENDS)}.")
        self._backend = backend_cls(**kwargs)

    def create_page(self, path, title, content, **kwargs):
        """See `base.WikiBackend.create_page` for the shared contract, and
        the chosen backend's own `create_page` docstring for what its
        backend-specific `**kwargs` accept."""
        return self._backend.create_page(path, title, content, **kwargs)
