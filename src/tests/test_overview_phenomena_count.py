# tests/test_overview_phenomena_count.py

"""
UX.110: the phenomena figure on the Home page (the overview) is the number of rows the Phenomena list shows, for a
galaxy of built phenomena and after every scatter, fill and edit; with the page cache on and off.
"""

import re

import pytest

from planetgen.db import countcache, dataversion
from planetgen.db import query as queryDb
from planetgen.db import store
from planetgen.generation.config import SystemConfig
from planetgen.generation.phenomena.nebula import Nebula
from planetgen.generation.phenomena.rogue import RoguePlanet
from planetgen.api.config import Config
from planetgen.web.app import create_app


@pytest.fixture(params=["cache-off", "cache-on"])
def site(request, mysql_config, monkeypatch):
    class RealConfig(Config):
        MYSQL_CONFIG = mysql_config
        WRITE_MYSQL_CONFIG = mysql_config
        CONTROL_MYSQL_CONFIG = mysql_config
        WEB_DATABASE = ""
        SESSION_COOKIE_SECURE = False
        SECRET_KEY = "test-secret"

    if request.param == "cache-off":
        monkeypatch.setenv("PLANETGEN_PAGE_CACHE", "off")
    application = create_app(RealConfig)
    application.testing = True
    return application.test_client()


def _save(mysql_config, phenomenon, kind, name):
    phenomenon.name = name
    saved = store.save_phenomenon(phenomenon, SystemConfig(), kind, config=mysql_config)
    dataversion.record(mysql_config.database, "generated")   # what a generation run does when it finishes
    return saved


def _overview_count(site):
    page = site.get("/").get_data(as_text=True)
    return int(re.search(r"<strong>([\d,]+)</strong> phenomen(?:on|a)", page).group(1).replace(",", ""))


def _list_total(site):
    page = site.get("/phenomena").get_data(as_text=True)
    return int(re.search(r'data-total="(\d+)"', page).group(1))


def _rows(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        return len(queryDb.list_phenomena(conn))
    finally:
        conn.close()


def _agree(site, mysql_config, expected):
    assert _rows(mysql_config) == expected
    assert _list_total(site) == expected
    assert _overview_count(site) == expected


def test_the_overview_count_follows_every_scatter_fill_and_edit(site, mysql_config):
    cfg = SystemConfig()
    _agree(site, mysql_config, 0)
    _save(mysql_config, Nebula(cfg), "nebula", "Alpha Cloud")
    _save(mysql_config, RoguePlanet(cfg), "rogue-planet", "Beta Wanderer")
    _agree(site, mysql_config, 2)

    # A scatter places five neutron stars no sector has built.
    conn = store.get_connection(mysql_config)
    rows = [(0, 0, slot, "neutron-star", None, slot, 0, 0, None, None, None, slot + 1, None, None) for slot in range(5)]
    store.insert_phenomenon_scatter(conn, rows)
    store.record_phenomenon_scatter_classes(conn, {("neutron-star", ""): 5})
    conn.commit()
    dataversion.record(mysql_config.database, "generated")
    _agree(site, mysql_config, 7)

    # A sector fill builds one of them: it leaves the unbuilt count and (as a star's own remnant, here a
    # standalone neutron star row) joins the built ones.
    first = conn.execute("SELECT MIN(id) AS i FROM phenomenon_scatter").fetchone()["i"]
    store.mark_phenomena_built(conn, [first])
    conn.commit()
    dataversion.record(mysql_config.database, "generated")
    _agree(site, mysql_config, 6)

    # An edit: a built phenomenon is renamed, and another is added.
    _save(mysql_config, Nebula(cfg), "nebula", "Zeta Cloud")
    _agree(site, mysql_config, 7)
    conn.close()


def test_a_capped_count_is_shown_capped_on_the_overview_as_in_the_list(site, mysql_config, monkeypatch, redis_server):
    """With no stored count yet the API's total stops at its cap (`total_capped`): the Phenomena list shows
    "9,999+", and the overview used to show the bare "9,999"."""
    cfg = SystemConfig()
    for n in range(4):
        _save(mysql_config, Nebula(cfg), "nebula", f"Cloud {n}")
    monkeypatch.setenv(countcache.ENV_VAR, "on")
    monkeypatch.setattr(countcache, "WAIT_SECONDS", 0.0)
    monkeypatch.setattr(countcache, "_start", lambda *a: None)   # no background count: the estimate answers
    monkeypatch.setattr(queryDb, "COUNT_FALLBACK_CAP", 3)
    countcache.forget_all()
    try:
        overview = site.get("/").get_data(as_text=True)
        shown = re.search(r"<strong>([\d,]+\+?)</strong> phenomen(?:on|a)", overview).group(1)
        listed = re.search(r'datatable-count-line"[^>]*>([\d,]+\+?) ', site.get("/phenomena").get_data(as_text=True)).group(1)
    finally:
        countcache.forget_all()
    assert shown == listed == "2+"
