# tests/test_wikijs_client.py

"""
Unit tests for `wikijs.client.WikiJsClient` -- every request/response cycle
is mocked at `urllib.request.urlopen`, so this file needs no live Wiki.js
instance and no network access to run. See `test_wikijs_client_integration.py`
for the optional live-instance counterpart.
"""

import json
import urllib.error
from unittest.mock import MagicMock, patch

import pytest

from wikijs import WikiJsAuthError, WikiJsClient, WikiJsPageExistsError, WikiJsRequestError

BASE_URL = "https://wiki.example.com"
API_TOKEN = "test-token-123"


def _fake_response(body_dict, status=200):
    """A context-manager-compatible stand-in for `urlopen`'s return value."""
    body_bytes = json.dumps(body_dict).encode("utf-8")
    response = MagicMock()
    response.__enter__ = MagicMock(return_value=response)
    response.__exit__ = MagicMock(return_value=False)
    response.read = MagicMock(return_value=body_bytes)
    response.status = status
    return response


def _create_result(succeeded=True, error_code=0, message="", page=None):
    return {
        "data": {
            "pages": {
                "create": {
                    "responseResult": {
                        "succeeded": succeeded,
                        "errorCode": error_code,
                        "slug": "",
                        "message": message,
                    },
                    "page": page,
                }
            }
        }
    }


@pytest.fixture
def client():
    return WikiJsClient(BASE_URL, API_TOKEN)


def test_create_page_success_sends_expected_request(client):
    page = {"id": 42, "path": "systems/kepler-442", "title": "Kepler-442", "updatedAt": "2026-01-01T00:00:00Z"}
    with patch("urllib.request.urlopen", return_value=_fake_response(_create_result(page=page))) as mock_urlopen:
        result = client.create_page(
            path="systems/kepler-442",
            title="Kepler-442",
            content="# Kepler-442\n\nA system.",
            description="A generated star system.",
            tags=["system", "generated"],
        )

    assert result.id == 42
    assert result.path == "systems/kepler-442"
    assert result.title == "Kepler-442"
    assert result.url == "https://wiki.example.com/systems/kepler-442"

    request = mock_urlopen.call_args[0][0]
    assert request.full_url == "https://wiki.example.com/graphql"
    assert request.get_header("Authorization") == f"Bearer {API_TOKEN}"
    assert request.get_header("Content-type") == "application/json"

    sent_body = json.loads(request.data.decode("utf-8"))
    variables = sent_body["variables"]
    assert variables["path"] == "systems/kepler-442"
    assert variables["title"] == "Kepler-442"
    assert variables["content"] == "# Kepler-442\n\nA system."
    assert variables["description"] == "A generated star system."
    assert variables["editor"] == "markdown"
    assert variables["locale"] == "en"
    assert variables["tags"] == ["system", "generated"]
    assert variables["isPublished"] is True
    assert variables["isPrivate"] is False
    assert "pages" in sent_body["query"]
    assert "create" in sent_body["query"]


def test_create_page_defaults_tags_to_empty_list(client):
    page = {"id": 1, "path": "p", "title": "t", "updatedAt": ""}
    with patch("urllib.request.urlopen", return_value=_fake_response(_create_result(page=page))) as mock_urlopen:
        client.create_page(path="p", title="t", content="c")

    request = mock_urlopen.call_args[0][0]
    variables = json.loads(request.data.decode("utf-8"))["variables"]
    assert variables["tags"] == []


def test_create_page_duplicate_path_raises_page_exists_error(client):
    result = _create_result(succeeded=False, error_code=2001, message="A page already exists at this path.")
    with patch("urllib.request.urlopen", return_value=_fake_response(result)):
        with pytest.raises(WikiJsPageExistsError):
            client.create_page(path="systems/kepler-442", title="Kepler-442", content="...")


def test_create_page_other_logical_failure_raises_request_error(client):
    result = _create_result(succeeded=False, error_code=9999, message="Something else went wrong.")
    with patch("urllib.request.urlopen", return_value=_fake_response(result)):
        with pytest.raises(WikiJsRequestError, match="Something else went wrong."):
            client.create_page(path="p", title="t", content="c")


def test_create_page_http_401_raises_auth_error(client):
    http_error = urllib.error.HTTPError(
        url=f"{BASE_URL}/graphql", code=401, msg="Unauthorized", hdrs=None, fp=None,
    )
    http_error.read = MagicMock(return_value=b'{"message": "invalid token"}')
    with patch("urllib.request.urlopen", side_effect=http_error):
        with pytest.raises(WikiJsAuthError):
            client.create_page(path="p", title="t", content="c")


def test_create_page_http_403_raises_auth_error(client):
    http_error = urllib.error.HTTPError(
        url=f"{BASE_URL}/graphql", code=403, msg="Forbidden", hdrs=None, fp=None,
    )
    http_error.read = MagicMock(return_value=b"")
    with patch("urllib.request.urlopen", side_effect=http_error):
        with pytest.raises(WikiJsAuthError):
            client.create_page(path="p", title="t", content="c")


def test_create_page_other_http_error_raises_request_error(client):
    http_error = urllib.error.HTTPError(
        url=f"{BASE_URL}/graphql", code=500, msg="Internal Server Error", hdrs=None, fp=None,
    )
    http_error.read = MagicMock(return_value=b"boom")
    with patch("urllib.request.urlopen", side_effect=http_error):
        with pytest.raises(WikiJsRequestError):
            client.create_page(path="p", title="t", content="c")


def test_create_page_unreachable_host_raises_request_error(client):
    with patch("urllib.request.urlopen", side_effect=urllib.error.URLError("nodename nor servname provided")):
        with pytest.raises(WikiJsRequestError, match="Could not reach Wiki.js"):
            client.create_page(path="p", title="t", content="c")


def test_create_page_malformed_json_raises_request_error(client):
    response = MagicMock()
    response.__enter__ = MagicMock(return_value=response)
    response.__exit__ = MagicMock(return_value=False)
    response.read = MagicMock(return_value=b"not json")
    with patch("urllib.request.urlopen", return_value=response):
        with pytest.raises(WikiJsRequestError, match="unparseable"):
            client.create_page(path="p", title="t", content="c")


def test_create_page_graphql_auth_error_raises_auth_error(client):
    body = {"errors": [{"message": "Access Denied"}]}
    with patch("urllib.request.urlopen", return_value=_fake_response(body)):
        with pytest.raises(WikiJsAuthError):
            client.create_page(path="p", title="t", content="c")


def test_create_page_graphql_other_error_raises_request_error(client):
    body = {"errors": [{"message": "Cannot query field 'bogus' on type 'Query'."}]}
    with patch("urllib.request.urlopen", return_value=_fake_response(body)):
        with pytest.raises(WikiJsRequestError):
            client.create_page(path="p", title="t", content="c")


def test_create_page_unexpected_shape_raises_request_error(client):
    with patch("urllib.request.urlopen", return_value=_fake_response({"data": {"pages": {}}})):
        with pytest.raises(WikiJsRequestError, match="Unexpected response shape"):
            client.create_page(path="p", title="t", content="c")
