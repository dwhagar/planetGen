# tests/test_wikiclient_client.py

"""
Unit tests for `wikiClient.client.WikiClient` -- the facade that dispatches
`create_page` to whichever backend it was constructed with. Each backend's
own request/response behavior is covered by `test_wikiclient_wikijs.py`/
`test_wikiclient_mediawiki.py`; this file only checks that `WikiClient`
picks the right backend class and forwards arguments/return values
untouched.
"""

from unittest.mock import MagicMock, patch

import pytest

from wikiClient import WikiClient
from wikiClient.base import WikiPage


def test_wikijs_backend_is_constructed_and_used():
    fake_page = WikiPage(id=1, path="p", title="t", url="https://wiki.example.com/p")
    mock_backend_cls = MagicMock()
    mock_backend_cls.return_value.create_page.return_value = fake_page

    # `WikiClient.__init__` looks the backend class up from `_BACKENDS`
    # (built once at module import time), so the mock has to replace that
    # dict entry directly -- patching the `WikiJsBackend` name in
    # `wikiClient.client`'s namespace wouldn't reach a reference `_BACKENDS`
    # already holds.
    with patch.dict("wikiClient.client._BACKENDS", {"wikijs": mock_backend_cls}):
        conduit = WikiClient(backend="wikijs", base_url="https://wiki.example.com", api_token="tok")
        result = conduit.create_page(path="p", title="t", content="c", tags=["x"])

    mock_backend_cls.assert_called_once_with(base_url="https://wiki.example.com", api_token="tok")
    mock_backend_cls.return_value.create_page.assert_called_once_with("p", "t", "c", tags=["x"])
    assert result is fake_page


def test_mediawiki_backend_is_constructed_and_used():
    fake_page = WikiPage(id=2, path="P", title="P", url="https://wiki.example.com/w/index.php?title=P")
    mock_backend_cls = MagicMock()
    mock_backend_cls.return_value.create_page.return_value = fake_page

    with patch.dict("wikiClient.client._BACKENDS", {"mediawiki": mock_backend_cls}):
        conduit = WikiClient(
            backend="mediawiki", base_url="https://wiki.example.com/w", username="Bot@x", password="pw"
        )
        result = conduit.create_page(path="P", title="P", content="c")

    mock_backend_cls.assert_called_once_with(base_url="https://wiki.example.com/w", username="Bot@x", password="pw")
    mock_backend_cls.return_value.create_page.assert_called_once_with("P", "P", "c")
    assert result is fake_page


def test_unknown_backend_raises_value_error():
    with pytest.raises(ValueError, match="Unknown wiki backend"):
        WikiClient(backend="dokuwiki", base_url="https://wiki.example.com")
