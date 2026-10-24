# tests/test_reset_counts.py

"""
UX.161: after `planetgen reset` the Home page's counts read the wiped database at once, with the stored counts
(`countcache`) and the page cache (`pagecache`) both on, as on a live site, and with them off.
"""

import re

import pytest

from planetgen.api.config import Config
from planetgen.cli.reset import reset_database
from planetgen.db import countcache
from planetgen.web.app import create_app
from tests.test_api import _place_sector


@pytest.fixture(params=["caches-off", "caches-on"])
def site(request, mysql_config, monkeypatch, redis_server):
    class RealConfig(Config):
        MYSQL_CONFIG = mysql_config
        WRITE_MYSQL_CONFIG = mysql_config
        CONTROL_MYSQL_CONFIG = mysql_config
        WEB_DATABASE = ""
        SESSION_COOKIE_SECURE = False
        SECRET_KEY = "test-secret"

    if request.param == "caches-off":
        monkeypatch.setenv("PLANETGEN_PAGE_CACHE", "off")
    else:
        monkeypatch.setenv(countcache.ENV_VAR, "on")
    countcache.forget_all()
    application = create_app(RealConfig)
    application.testing = True
    yield application.test_client()
    countcache.forget_all()


def _home_counts(site):
    page = site.get("/").get_data(as_text=True)
    return {noun: int(number.replace(",", ""))
            for number, noun in re.findall(r"<strong>([\d,]+)\+?</strong> ([a-z ]+?)</a>", page)}


def test_the_home_counts_are_zero_right_after_a_reset(site, mysql_config):
    _place_sector(mysql_config, "One", (5.0, 5.0, 5.0))
    _place_sector(mysql_config, "Two", (50.0, 5.0, 5.0))
    before = _home_counts(site)
    assert before["sectors"] == 2 and before["systems"] == 2
    assert _home_counts(site) == before   # served again: the page cache and the stored counts hold it

    assert reset_database(mysql_config, assume_yes=True)
    after = _home_counts(site)
    assert after and not any(after.values()), after
