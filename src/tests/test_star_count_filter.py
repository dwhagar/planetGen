# tests/test_star_count_filter.py

"""UX.239: the stars filter of the search (exactly N or N or more; Single, Binary, Trinary, Quaternary, Quintenary and beyond) in `util/star_count`, `db.query.search`,
`GET /api/search?stars=` and the search page and its box (`stars:3`)."""

import pytest

from planetgen.db import store
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.util import draw, star_count
from tests.test_api import client, seeded_sector  # noqa: F401


# --- the parser ---------------------------------------------------------------------------------------------------------------

@pytest.mark.parametrize("text,expected", [
    ("", None), ("  ", None), (None, None),
    ("1", (1, 1)), ("2", (2, 2)), ("3", (3, 3)), ("12", (12, 12)),
    ("3+", (3, None)), ("1+", (1, None)), (" 4 + ", (4, None)),
    ("single", (1, 1)), ("binary", (2, 2)), ("trinary", (3, 3)), ("quaternary", (4, 4)),
    ("quintenary", (5, None)), ("Quintenary+", (5, None)), ("BINARY", (2, 2)), ("binary+", (2, None)),
])
def test_the_stars_filter_reads_counts_names_and_or_more(text, expected):
    assert star_count.parse_star_count(text) == expected


@pytest.mark.parametrize("text", ["0", "0+", "-1", "three", "3-4", "1.5", "stars", "+", "3++"])
def test_a_stars_filter_that_is_none_of_these_is_refused(text):
    with pytest.raises(ValueError, match="stars must be"):
        star_count.parse_star_count(text)


def test_every_choice_the_form_offers_parses_and_every_name_is_a_choice_or_the_count_of_one():
    for value, label in star_count.STAR_COUNT_CHOICES:
        assert star_count.parse_star_count(value) is not None and label
    assert [label for _v, label in star_count.STAR_COUNT_CHOICES] == ["Single", "Binary", "Trinary", "Quaternary", "Quintenary and beyond"]
    assert {name for name, _n in star_count.STAR_COUNT_NAMES} == {"single", "binary", "trinary", "quaternary", "quintenary"}


def test_the_labels_say_the_choice_in_words():
    assert star_count.star_count_label((2, 2)) == "Binary"
    assert star_count.star_count_label((5, None)) == "Quintenary and beyond"
    assert star_count.star_count_label((7, None)) == "7 stars or more"
    assert star_count.star_count_label((6, 6)) == "6 stars"
    assert star_count.star_count_label(None) == ""


# --- the search --------------------------------------------------------------------------------------------------------------

def _system(mysql_config, name, stars):
    config = SystemConfig()
    config.PLANETS = False
    config.STARS = stars
    with draw.bound(f"star-count-filter:{name}"):
        system = StarSystem(system_config=config)
    assert len(system.stars) == stars, (name, len(system.stars))
    system.name = name
    return store.save_system(system, config, config=mysql_config)


@pytest.fixture
def multiples(client, seeded_sector):  # noqa: F811
    mysql_config, _sector_id, _systems = seeded_sector
    for name, stars in (("Solo Kappa", 1), ("Pair Kappa", 2), ("Triple Kappa", 3), ("Quad Kappa", 4)):
        _system(mysql_config, name, stars)
    return client


def _found(response):
    assert response.status_code == 200, response.get_data(as_text=True)
    return {row["name"]: row for row in response.get_json()["results"]["systems"]["rows"]}


def test_the_stars_filter_finds_a_known_binary_and_a_known_triple(multiples):
    binary = _found(multiples.get("/api/search?stars=2&panels=systems"))
    triple = _found(multiples.get("/api/search?stars=trinary&panels=systems"))
    assert "Pair Kappa" in binary and binary["Pair Kappa"]["star_count"] == 2
    assert "Triple Kappa" in triple and triple["Triple Kappa"]["star_count"] == 3
    assert all(row["star_count"] == 2 for row in binary.values()) and all(row["star_count"] == 3 for row in triple.values())
    assert "Triple Kappa" not in binary and "Pair Kappa" not in triple


def test_or_more_and_the_named_choices(multiples):
    three_or_more = _found(multiples.get("/api/search?stars=3%2B&panels=systems"))
    assert {"Triple Kappa", "Quad Kappa"} <= set(three_or_more) and not {"Solo Kappa", "Pair Kappa"} & set(three_or_more)
    assert all(row["star_count"] >= 3 for row in three_or_more.values())
    assert "Solo Kappa" in _found(multiples.get("/api/search?stars=single&panels=systems"))
    assert "Quad Kappa" in _found(multiples.get("/api/search?stars=quaternary&panels=systems"))
    assert not {"Triple Kappa", "Quad Kappa"} & set(_found(multiples.get("/api/search?stars=quintenary&panels=systems")))


def test_the_filter_and_a_name_together_keep_only_what_matches_both(multiples):
    both = _found(multiples.get("/api/search?stars=3%2B&system_q=Kappa&panels=systems"))
    assert set(both) == {"Triple Kappa", "Quad Kappa"}
    assert not _found(multiples.get("/api/search?stars=2&system_q=Triple&panels=systems"))


def test_a_stars_filter_that_is_not_one_is_a_400(multiples):
    response = multiples.get("/api/search?stars=lots")
    assert response.status_code == 400 and "stars must be" in response.get_json()["error"]


def test_without_the_filter_the_systems_panel_still_needs_a_name(multiples):
    assert multiples.get("/api/search?panels=systems").get_json()["results"]["systems"] is None
    assert multiples.get("/api/search?stars=2&panels=none").get_json()["results"]["systems"] is None
