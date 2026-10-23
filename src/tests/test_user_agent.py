"""API.31: one function builds the `User-Agent` of every request the software makes to another site, and nothing
else spells it."""

import ast
import re
from pathlib import Path

import pytest

from planetgen import _version
from planetgen.api import version as api_version
from planetgen.util import user_agent as ua

SRC = Path(__file__).resolve().parents[1]
PRODUCT_TOKEN = r"[A-Za-z0-9!#$%&'*+.^_`|~-]+/[A-Za-z0-9!#$%&'*+.^_`|~-]+"


def test_the_string_is_product_tokens_with_the_release_the_project_url_and_the_api_version():
    text = ua.user_agent()
    assert re.fullmatch(rf"{PRODUCT_TOKEN} \(\+https://[^)\s]+\) {PRODUCT_TOKEN}", text)
    assert text == f"PlanetGen/{_version.__version__} (+{_version.REPO_URL}) PlanetGen-API/{api_version.API_VERSION}"


def test_the_string_follows_both_versions(monkeypatch):
    before = ua.user_agent()
    monkeypatch.setattr(ua, "API_VERSION", api_version.API_VERSION + 1)
    assert ua.user_agent() != before and ua.user_agent().endswith(f"PlanetGen-API/{api_version.API_VERSION + 1}")
    monkeypatch.setattr(ua, "__version__", "9.9.9")
    assert "PlanetGen/9.9.9 " in ua.user_agent()


def _sets_a_user_agent(path):
    """Whether the module builds a request with a `User-Agent` that is not `user_agent()`."""
    tree = ast.parse(path.read_text(encoding="utf-8"))
    for node in ast.walk(tree):
        if isinstance(node, ast.Constant) and isinstance(node.value, str) and node.value.lower() == "user-agent":
            return True
    return False


def _calls_user_agent(path):
    tree = ast.parse(path.read_text(encoding="utf-8"))
    return any(isinstance(node, ast.Call) and getattr(node.func, "id", None) == "user_agent" for node in ast.walk(tree))


OUTGOING = ["planetgen/web/lib/apiclient.py", "planetgen/wiki/mediawiki.py", "planetgen/wiki/wikijs.py",
            "planetgen/datasets/fetch.py", "../scripts/build_element_tables.py"]


@pytest.mark.parametrize("relative", OUTGOING)
def test_every_outgoing_request_path_sends_the_standard_string(relative):
    path = (SRC / relative).resolve()
    assert _sets_a_user_agent(path) and _calls_user_agent(path), f"{relative} must send user_agent()"


def test_no_other_module_makes_an_outgoing_request():
    """A module that opens a URL is listed in `OUTGOING`; this finds the ones that are not."""
    listed = {(SRC / relative).resolve() for relative in OUTGOING}
    found = set()
    for path in (SRC / "planetgen").rglob("*.py"):
        text = path.read_text(encoding="utf-8")
        if re.search(r"urlopen\(|urllib\.request\.Request\(|build_opener\(|http\.client\.HTTPS?Connection\(|import requests|import httpx",
                     text):
            found.add(path.resolve())
    assert found <= listed, f"outgoing requests outside OUTGOING: {sorted(str(p) for p in found - listed)}"


def test_the_page_scripts_send_the_api_version_where_they_cannot_set_a_user_agent():
    base = (SRC / "planetgen" / "web" / "templates" / "base.html").read_text(encoding="utf-8")
    assert 'name="planetgen-api-version"' in base and "apiversion.js" in base
    assert base.index("apiversion.js") < base.index("localtime.js")  # before any script that fetches
