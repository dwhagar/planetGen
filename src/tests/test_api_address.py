# tests/test_api_address.py

"""
OPS.86: the API address is the one in config.json (`api_base_url`) and nothing falls back to a loopback address: the
setting has no default, an empty one is refused by the HTTP transport with a plain message, and the example config
does not carry a loopback address as if it were the answer.
"""

import json
import pathlib

import pytest

from planetgen.util import settings as settings_module
from planetgen.web.lib import apiclient

ROOT = pathlib.Path(__file__).resolve().parents[2]


def test_api_base_url_has_no_default_and_no_loopback_address():
    assert settings_module.Settings.model_fields["api_base_url"].default == ""
    schema = json.loads((ROOT / "config.schema.json").read_text())
    assert schema["properties"]["api_base_url"]["default"] == ""
    assert "127.0.0.1" not in (ROOT / "config.json.example").read_text().split('"api_base_url"', 1)[1].split("\n", 1)[0]


@pytest.mark.parametrize("value", ["https://example.org/api", "http://127.0.0.1:8000/api", ""])
def test_api_base_url_accepts_an_address_or_nothing(value):
    assert settings_module.Settings(api_base_url=value).api_base_url == value


def test_api_base_url_drops_a_trailing_slash_and_refuses_what_is_not_an_address():
    assert settings_module.Settings(api_base_url="https://example.org/api/").api_base_url == "https://example.org/api"
    with pytest.raises(ValueError):
        settings_module.Settings(api_base_url="example.org/api")


def test_the_http_transport_refuses_with_a_plain_message_when_no_address_is_set(monkeypatch):
    monkeypatch.setattr(apiclient, "API_BASE_URL", "")
    monkeypatch.setattr(apiclient, "_transport", None)

    def never(*args, **kwargs):
        raise AssertionError("a connection was attempted with no address configured")

    monkeypatch.setattr(apiclient.urllib.request, "urlopen", never)
    with pytest.raises(apiclient.ApiError, match="No API address is configured: set api_base_url in config.json"):
        apiclient._http_transport("GET", "/health", None, {}, 5)
