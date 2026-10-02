# tests/test_search_edge_cases.py

"""
Full-text search edge cases (TEST.18): the name search behind `/search`
(`web/views.py`), `/api/search` (`api/routes.py`) and `queryDb.search`,
whose whole-word matching (`queryDb._name_match`) splits the typed text
into words for `MATCH ... AGAINST` in boolean mode, with a REGEXP for
the words the FULLTEXT index leaves out. Covered: boolean-mode operator
characters in names and terms, words shorter than the server's
`innodb_ft_min_token_size` or longer than its `innodb_ft_max_token_size`,
the server's own stopwords, and LIKE wildcards typed into a search.

Token sizes and stopwords are read from the server at run time, so each
engine (MySQL 8, MariaDB) is held to its own settings. Every test takes
`mysql_config` (see `conftest.py`).
"""

from urllib.parse import urlencode

import markupsafe
import pytest

import queryDb
from stellarObjects import _db
from stellarObjects.config import SystemConfig
from stellarObjects.galaxyGeometry import sector_position_pc
from stellarObjects.spaceSector import SpaceSector
from stellarObjects.systemData import StarSystem

from tests.test_web_pages import db_client  # noqa: F401 -- fixture
from tests.test_web_search import _panel


def _insert_sectors(config, names, addressed=False):
    conn = _db.get_connection(config)
    try:
        with conn:
            for slot, name in enumerate(names):
                if addressed:
                    x, y, z = sector_position_pc(5, 0, slot, 4.0)
                    conn.execute(
                        "INSERT INTO sectors (name, edge_mpc, center_x_pc, center_y_pc, center_z_pc,"
                        " galactic_radius_pc, ring_index, layer_index, ring_slot_index)"
                        " VALUES (?, 4000, ?, ?, ?, ?, 5, 0, ?)",
                        (name, x, y, z, (x * x + y * y) ** 0.5, slot),
                    )
                else:
                    conn.execute("INSERT INTO sectors (name, edge_mpc) VALUES (?, 4000)", (name,))
    finally:
        conn.close()


def _server_value(config, sql):
    conn = _db.get_connection(config)
    try:
        return conn.execute(sql).fetchone()["v"]
    finally:
        conn.close()


def _sector_names(config, term):
    conn = queryDb.open_readonly(config)
    try:
        texts = {"sector_q": term}
        result = queryDb.search(conn, texts, {facet: set() for facet in queryDb.SEARCH_TAG_FACETS})
        return sorted(row["name"] for row in result["results"]["sectors"]["rows"])
    finally:
        conn.close()


def _server_stopwords(config):
    """The stopwords the server's InnoDB FULLTEXT indexes leave out: its
    configured stopword table, else the built-in list."""
    conn = _db.get_connection(config)
    try:
        if not conn.execute("SELECT @@innodb_ft_enable_stopword AS v").fetchone()["v"]:
            return set()
        table = conn.execute("SELECT @@innodb_ft_server_stopword_table AS v").fetchone()["v"]
        if table:
            db, name = table.split("/", 1)
            rows = conn.execute(f"SELECT value FROM `{db}`.`{name}`").fetchall()
        else:
            rows = conn.execute("SELECT value FROM INFORMATION_SCHEMA.INNODB_FT_DEFAULT_STOPWORD").fetchall()
        return {row["value"].lower() for row in rows}
    finally:
        conn.close()


# --- Operator characters ------------------------------------------------------------

_OPERATOR_NAMES = [
    "O'Brien Reach", "Brien Hollow", "Kepler-42", "Kepler 420", 'Tarn "Quoted" Deep', "Vega* Prime",
    "Vega Outpost", "Vegatron", "Paren (Inner) Gap", "Angle <Less> Rim", "Tilde~Wave", "Mail@Host Drift",
    "Plus+Minus Gate",
]


@pytest.mark.parametrize("term, expected", [
    ("O'Brien", ["O'Brien Reach"]),
    ("brien", ["Brien Hollow", "O'Brien Reach"]),
    ("Kepler-42", ["Kepler-42"]),
    ('"Quoted"', ['Tarn "Quoted" Deep']),
    ('"tarn deep"', ['Tarn "Quoted" Deep']),
    # Typed operators are not boolean-mode syntax: -Prime still needs Prime.
    ("+Vega -Prime", ["Vega* Prime"]),
    ("Vega*", ["Vega Outpost", "Vega* Prime"]),
    ("(Inner)", ["Paren (Inner) Gap"]),
    ("<Less>", ["Angle <Less> Rim"]),
    ("~wave", ["Tilde~Wave"]),
    ("mail@host", ["Mail@Host Drift"]),
    ("Plus+Minus", ["Plus+Minus Gate"]),
    ("+-*\"'@()<>~", []),
    (">", []),
])
def test_operator_characters_in_names_and_terms(mysql_config, term, expected):
    _insert_sectors(mysql_config, _OPERATOR_NAMES)
    assert _sector_names(mysql_config, term) == expected


def test_operator_characters_through_the_api(db_client, mysql_config):
    _insert_sectors(mysql_config, _OPERATOR_NAMES)

    def names(**params):
        resp = db_client.get("/api/search?" + urlencode(params))
        assert resp.status_code == 200
        return sorted(row["name"] for row in resp.get_json()["results"]["sectors"]["rows"])

    assert names(sector_q="O'Brien") == ["O'Brien Reach"]
    assert names(sector_q='"Quoted" (Inner)') == []
    assert names(sector_q="+Vega -Prime") == ["Vega* Prime"]
    assert names(sector_q="<Less> ~") == ["Angle <Less> Rim"]
    assert names(sector_q="@") == []


def test_operator_characters_on_every_object_panel(db_client, mysql_config):
    sector = SpaceSector("Quiet Sector", edge_ly=11.5)
    cfg = SystemConfig()
    cfg.STAR_TYPE = "K1V"
    cfg.BINARY_SYSTEM = False
    cfg.PLANETS = True
    cfg.MOONS = True
    for _attempt in range(20):
        system = StarSystem(system_config=cfg)
        if any(getattr(planet, "moons", None) for planet in system.planets):  # belts have no moons
            break
    sector.add_system(system, position=(0.0, 0.0, 0.0), system_config=cfg)
    _db.save_sector(sector, config=mysql_config)
    names = {
        "star_systems": "Tarn's (Hope)", "stars": "Tarn's Hope A+", "planets": "Tarn's Hope-7 <b>",
        "moons": '"Ash" ~Minor',
    }
    conn = _db.get_connection(mysql_config)
    try:
        with conn:
            for table, name in names.items():
                conn.execute(f"UPDATE {table} SET name = ? ORDER BY id LIMIT 1", (name,))
    finally:
        conn.close()

    html = db_client.get("/search?" + urlencode({"q": "tarn's hope"})).get_data(as_text=True)
    for panel, table in (("systems", "star_systems"), ("stars", "stars"), ("planets", "planets")):
        assert f">{markupsafe.escape(names[table])}<" in _panel(html, panel)
    assert "search-moons" not in html or "Ash" not in _panel(html, "moons")

    html = db_client.get("/search?" + urlencode({"moon_q": '"ash" minor'})).get_data(as_text=True)
    moons = _panel(html, "moons")
    assert "(1)" in moons and str(markupsafe.escape(names["moons"])) in moons

    html = db_client.get("/search?" + urlencode({"planet_q": "Hope-7 <b>"})).get_data(as_text=True)
    assert "(1)" in _panel(html, "planets")
    assert "<b>" not in _panel(html, "planets")  # escaped, not markup


# --- Token sizes ------------------------------------------------------------------------

def test_words_shorter_than_the_min_token_size(mysql_config):
    min_size = _server_value(mysql_config, "SELECT @@innodb_ft_min_token_size AS v")
    if min_size < 2:
        pytest.skip(f"innodb_ft_min_token_size is {min_size}: every word is in the index")
    short = "Zqxw"[:min_size - 1]
    exact = "Korvax"[:min_size]
    _insert_sectors(mysql_config, [
        f"Ossiran {short}", f"{short}a Ossiran", "Ossiran", f"{exact} Station", f"{exact}x Station",
        "Tobar B", "Tobar Beta",
    ])
    assert _sector_names(mysql_config, short) == [f"Ossiran {short}"]
    assert _sector_names(mysql_config, short.lower() + " ossiran") == [f"Ossiran {short}"]
    assert _sector_names(mysql_config, exact) == [f"{exact} Station"]
    assert _sector_names(mysql_config, "b") == ["Tobar B"]
    assert _sector_names(mysql_config, "tobar b") == ["Tobar B"]


def test_words_longer_than_the_max_token_size(mysql_config):
    # A name can be 255 characters; a word past innodb_ft_max_token_size
    # (84 by default) isn't in the index, so MATCH alone never finds it.
    max_size = _server_value(mysql_config, "SELECT @@innodb_ft_max_token_size AS v")
    word = ("Llanfair" * 40)[:max_size + 1]
    _insert_sectors(mysql_config, [f"{word} Reach", f"{word}x Reach", "Plain Reach"])
    assert _sector_names(mysql_config, word) == [f"{word} Reach"]
    assert _sector_names(mysql_config, f"reach {word.lower()}") == [f"{word} Reach"]
    assert _sector_names(mysql_config, "reach") == sorted([f"{word} Reach", f"{word}x Reach", "Plain Reach"])


# --- Stopwords --------------------------------------------------------------------------

def test_every_server_stopword_is_searchable(mysql_config):
    stopwords = sorted(_server_stopwords(mysql_config))
    if not stopwords:
        pytest.skip("The server's FULLTEXT indexes keep every word")
    _insert_sectors(mysql_config, [f"Velorn {word.title()}" for word in stopwords] + ["Velorn"])
    missing = [
        word for word in stopwords
        if _sector_names(mysql_config, word) != [f"Velorn {word.title()}"]
        or _sector_names(mysql_config, f"velorn {word}") != [f"Velorn {word.title()}"]
    ]
    assert not missing, f"Stopwords this server's index leaves out that search can't find: {missing}"


# --- LIKE wildcards -----------------------------------------------------------------------

_WILDCARD_NAMES = ["100% Reach", "Plain Reach", "Under_Score Gate", "UnderXScore Gate", "Back\\Slash Rim",
                   "Backslash Rim"]


@pytest.mark.parametrize("term, expected", [
    ("%", []),
    ("_", []),
    ("\\", []),
    ("%%", []),
    ("100%", ["100% Reach"]),
    ("%reach%", ["100% Reach", "Plain Reach"]),
    ("Under_Score", ["Under_Score Gate"]),
    ("Under%Score", []),
    ("Back\\Slash", ["Back\\Slash Rim"]),
    ("Back\\", ["Back\\Slash Rim"]),  # "Back" as a whole word
])
def test_like_wildcards_typed_into_search(db_client, mysql_config, term, expected):
    _insert_sectors(mysql_config, _WILDCARD_NAMES)
    resp = db_client.get("/search?" + urlencode({"sector_q": term}))
    assert resp.status_code == 200
    html = resp.get_data(as_text=True)
    if expected:
        panel = _panel(html, "sectors")
        assert f"({len(expected)})" in panel
        found = [name for name in _WILDCARD_NAMES if f">{markupsafe.escape(name)}<" in panel]
        assert found == [name for name in _WILDCARD_NAMES if name in expected]
    else:
        assert not any(f">{markupsafe.escape(name)}<" in html for name in _WILDCARD_NAMES)


def test_like_wildcards_typed_into_the_galaxy_map_locate(db_client, mysql_config):
    # Its single short word matches the start of a name with LIKE.
    _insert_sectors(mysql_config, ["U_Gate", "UxGate", "B\\Rim", "BxRim", "B%Rim", "C Rim"], addressed=True)

    def names(term):
        resp = db_client.get("/api/galaxy/locate?" + urlencode({"q": term}))
        assert resp.status_code == 200
        return sorted(match["name"] for match in resp.get_json()["matches"])

    assert names("U_") == ["U_Gate"]
    assert names("B\\") == ["B\\Rim"]
    assert names("B%") == ["B%Rim"]
    assert names("B") == ["B%Rim", "B\\Rim", "BxRim"]
    assert names("%") == [] and names("_") == [] and names("\\") == []
