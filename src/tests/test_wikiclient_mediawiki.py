# tests/test_wikiclient_mediawiki.py

"""
Unit tests for `wikiClient.mediawiki.MediaWikiBackend` -- every request/
response cycle is mocked at the backend's own cookie-jar opener (`opener.
open`), so this file needs no live MediaWiki instance and no network access
to run. See `test_wikiclient_mediawiki_integration.py` for the optional
live-instance counterpart.
"""

import json
import urllib.error
from unittest.mock import MagicMock, patch

import pytest

from wikiClient import WikiClientAuthError, WikiClientPageExistsError, WikiClientRequestError
from wikiClient.mediawiki import MediaWikiBackend

BASE_URL = "https://wiki.example.com/w"
USERNAME = "Bot@planetgen"
PASSWORD = "test-password-123"

_LOGIN_TOKEN_RESPONSE = {"query": {"tokens": {"logintoken": "logintok+\\"}}}
_LOGIN_SUCCESS_RESPONSE = {"login": {"result": "Success", "lguserid": 1, "lgusername": "Bot"}}
_LOGIN_FAILURE_RESPONSE = {"login": {"result": "Failed", "reason": "Incorrect password entered."}}
_CSRF_TOKEN_RESPONSE = {"query": {"tokens": {"csrftoken": "csrftok+\\"}}}


def _fake_response(body_dict):
    """A context-manager-compatible stand-in for `opener.open`'s return value."""
    body_bytes = json.dumps(body_dict).encode("utf-8")
    response = MagicMock()
    response.__enter__ = MagicMock(return_value=response)
    response.__exit__ = MagicMock(return_value=False)
    response.read = MagicMock(return_value=body_bytes)
    return response


def _edit_success_response(page_id=7, title="Systems/Kepler-442"):
    return {"edit": {"result": "Success", "pageid": page_id, "title": title, "newrevid": 99}}


def _login_then(*bodies):
    """The 4-call sequence every first `create_page` makes: login token,
    login, csrf token, then whatever `bodies` gives for the edit call(s)
    after that."""
    return [_LOGIN_TOKEN_RESPONSE, _LOGIN_SUCCESS_RESPONSE, _CSRF_TOKEN_RESPONSE, *bodies]


@pytest.fixture
def client():
    return MediaWikiBackend(BASE_URL, USERNAME, PASSWORD)


def test_create_page_success_sends_expected_requests(client):
    responses = [_fake_response(body) for body in _login_then(_edit_success_response())]
    with patch.object(client._opener, "open", side_effect=responses) as mock_open:
        result = client.create_page(path="Systems/Kepler-442", title="Kepler-442", content="A system.")

    assert result.id == 7
    assert result.path == "Systems/Kepler-442"
    assert result.title == "Systems/Kepler-442"
    assert result.url.startswith(BASE_URL)

    assert mock_open.call_count == 4
    edit_request = mock_open.call_args_list[3][0][0]
    assert edit_request.full_url == f"{BASE_URL}/api.php"
    sent_params = dict(pair.split("=", 1) for pair in edit_request.data.decode("utf-8").split("&"))
    assert sent_params["title"] == "Systems%2FKepler-442"
    assert sent_params["createonly"] == "1"


def test_create_page_reuses_session_on_second_call(client):
    responses = [
        *[_fake_response(body) for body in _login_then(_edit_success_response())],
        _fake_response(_CSRF_TOKEN_RESPONSE),
        _fake_response(_edit_success_response(page_id=8, title="Systems/Other")),
    ]
    with patch.object(client._opener, "open", side_effect=responses) as mock_open:
        client.create_page(path="Systems/Kepler-442", title="Kepler-442", content="A system.")
        client.create_page(path="Systems/Other", title="Other", content="Another system.")

    # Second call skips login-token + login (already logged in), only
    # fetches a fresh csrf token and sends the edit -- 4 + 2 = 6 total.
    assert mock_open.call_count == 6


def test_create_page_login_failure_raises_auth_error(client):
    responses = [_fake_response(_LOGIN_TOKEN_RESPONSE), _fake_response(_LOGIN_FAILURE_RESPONSE)]
    with patch.object(client._opener, "open", side_effect=responses):
        with pytest.raises(WikiClientAuthError, match="Incorrect password"):
            client.create_page(path="p", title="t", content="c")


def test_create_page_article_exists_raises_page_exists_error(client):
    error_body = {"error": {"code": "articleexists", "info": "The article you tried to create has been created already"}}
    responses = [_fake_response(body) for body in _login_then(error_body)]
    with patch.object(client._opener, "open", side_effect=responses):
        with pytest.raises(WikiClientPageExistsError):
            client.create_page(path="p", title="t", content="c")


def test_create_page_permission_denied_raises_auth_error(client):
    error_body = {"error": {"code": "permissiondenied", "info": "You do not have permission to edit this page."}}
    responses = [_fake_response(body) for body in _login_then(error_body)]
    with patch.object(client._opener, "open", side_effect=responses):
        with pytest.raises(WikiClientAuthError):
            client.create_page(path="p", title="t", content="c")


def test_create_page_other_error_raises_request_error(client):
    error_body = {"error": {"code": "badtoken", "info": "Invalid token."}}
    responses = [_fake_response(body) for body in _login_then(error_body)]
    with patch.object(client._opener, "open", side_effect=responses):
        with pytest.raises(WikiClientRequestError, match="badtoken"):
            client.create_page(path="p", title="t", content="c")


def test_create_page_http_401_raises_auth_error(client):
    http_error = urllib.error.HTTPError(
        url=f"{BASE_URL}/api.php", code=401, msg="Unauthorized", hdrs=None, fp=None,
    )
    http_error.read = MagicMock(return_value=b"")
    with patch.object(client._opener, "open", side_effect=http_error):
        with pytest.raises(WikiClientAuthError):
            client.create_page(path="p", title="t", content="c")


def test_create_page_unreachable_host_raises_request_error(client):
    with patch.object(client._opener, "open", side_effect=urllib.error.URLError("nodename nor servname provided")):
        with pytest.raises(WikiClientRequestError, match="Could not reach MediaWiki"):
            client.create_page(path="p", title="t", content="c")


def test_create_page_malformed_json_raises_request_error(client):
    response = MagicMock()
    response.__enter__ = MagicMock(return_value=response)
    response.__exit__ = MagicMock(return_value=False)
    response.read = MagicMock(return_value=b"not json")
    with patch.object(client._opener, "open", return_value=response):
        with pytest.raises(WikiClientRequestError, match="unparseable"):
            client.create_page(path="p", title="t", content="c")


def test_create_page_unexpected_shape_raises_request_error(client):
    responses = [_fake_response(body) for body in _login_then({"edit": {"result": "Success"}})]
    with patch.object(client._opener, "open", side_effect=responses):
        with pytest.raises(WikiClientRequestError, match="Unexpected response shape"):
            client.create_page(path="p", title="t", content="c")
