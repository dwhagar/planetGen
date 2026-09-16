# wikijs/__init__.py

"""
A small, standalone GraphQL client for publishing pages to a Wiki.js
instance -- see `client.py`'s module docstring for the full design.

    from wikijs import WikiJsClient, WikiJsPageExistsError

    client = WikiJsClient(base_url="https://wiki.example.com", api_token="...")
    page = client.create_page(path="systems/kepler-442", title="Kepler-442", content="...")
"""

from .client import WikiJsClient, WikiPage
from .exceptions import WikiJsAuthError, WikiJsError, WikiJsPageExistsError, WikiJsRequestError

__all__ = [
    "WikiJsClient",
    "WikiPage",
    "WikiJsError",
    "WikiJsAuthError",
    "WikiJsPageExistsError",
    "WikiJsRequestError",
]
