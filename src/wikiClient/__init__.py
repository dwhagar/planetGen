# wikiClient/__init__.py

"""
A single library for publishing pages to either a Wiki.js or a MediaWiki
instance, behind one `WikiClient` object -- see `client.py`'s module
docstring for the full design, and `wikijs.py`/`mediawiki.py` for each
backend's own notes.

    from wikiClient import WikiClient, WikiClientPageExistsError

    client = WikiClient(backend="wikijs", base_url="https://wiki.example.com", api_token="...")
    page = client.create_page(path="systems/kepler-442", title="Kepler-442", content="...")
"""

from .base import WikiPage
from .client import WikiClient
from .exceptions import WikiClientAuthError, WikiClientError, WikiClientPageExistsError, WikiClientRequestError
from .mediawiki import MediaWikiBackend
from .wikijs import WikiJsBackend

__all__ = [
    "WikiClient",
    "WikiPage",
    "WikiClientError",
    "WikiClientAuthError",
    "WikiClientPageExistsError",
    "WikiClientRequestError",
    "WikiJsBackend",
    "MediaWikiBackend",
]
