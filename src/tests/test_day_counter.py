# tests/test_day_counter.py

"""UX.243: the day-and-night counter of an object page: what `web.lib.day_counter.counter` makes of a stored row, and the panel on the page of a planet
and a moon (the clock's own numbers are in `tests/js/daycounter.test.mjs`)."""

import re

import pytest

from planetgen.physics import day_night
from planetgen.web.lib import day_counter
from tests import test_item_pages as pages
from tests.test_item_pages import MOON, PLANET, STAR, client  # noqa: F401


RING_STATES = ("permanent day side", "permanent day side", "permanent day side", "terminator band", "permanent night side", "permanent night side",
               "permanent night side")


def _row(state=day_night.CYCLE, **extra):
    row = {"day_night_state": state, "solar_day_hours": 24.0, "night_hours": 12.0, "days_per_year": 365.25}
    for number in range(1, 8):
        row.update({f"band{number}_mean_k": 250.0 + number, f"band{number}_day_hours": 12.0, f"band{number}_night_hours": 12.0,
                    f"band{number}_state": "day and night"})
    row.update(extra)
    return row


def test_a_body_with_no_stored_state_has_no_counter():
    assert day_counter.counter({"solar_day_hours": 24.0}, "planet") is None
    assert day_counter.counter(None, "planet") is None


def test_a_body_with_a_day_ticks_and_carries_the_stored_hours():
    made = day_counter.counter(_row(), "planet")
    assert made["ticks"] and made["day_hours"] == 24.0 and made["night_hours"] == 12.0 and made["epoch_ms"] == day_counter.EPOCH_MS
    assert made["text"] == day_counter.STATE_TEXT[day_night.CYCLE]
    assert "faster than the real rotation" in made["note"]


def test_a_tidally_locked_body_states_its_permanent_sides_and_counts_no_day():
    for kind, toward in (("planet", "its star"), ("moon", "its planet")):
        made = day_counter.counter(_row(day_night.LOCKED), kind)
        assert not made["ticks"] and made["day_hours"] is None
        assert "Tidally locked" in made["text"] and toward in made["text"]
    assert "permanent day" in day_counter.counter(_row(day_night.LOCKED), "planet")["text"]


def test_every_state_word_has_a_sentence():
    for state in day_night.STATES:
        assert day_counter.counter(_row(state), "planet")["text"]


def test_a_body_with_no_solar_day_ticks_nothing():
    assert not day_counter.counter(_row(solar_day_hours=None), "planet")["ticks"]
    assert not day_counter.counter(_row(solar_day_hours=0.0), "planet")["ticks"]


def test_the_seven_bands_show_their_temperature_and_their_day_and_night_beside_it():
    made = day_counter.counter(_row(), "planet")
    assert [band["label"] for band in made["bands"]] == ["North polar", "North middle", "North tropics", "Equator", "South tropics", "South middle", "South polar"]
    first = made["bands"][0]
    assert 'data-q="temperature"' in str(first["cells"][0]) and 'data-q="time"' in str(first["cells"][1]) and 'data-q="time"' in str(first["cells"][2])
    assert made["band_headings"] == ["Annual mean temperature", "Day", "Night", "State"]


def test_a_world_that_turns_once_an_orbit_shows_its_rings():
    row = _row(day_night.LOCKED, **{f"band{n}_state": s for n, s in enumerate(RING_STATES, start=1)})
    assert [band["label"] for band in day_counter.counter(row, "planet")["bands"]][0] == "Day centre"


def test_a_body_with_no_bands_shows_none():
    assert day_counter.counter({"day_night_state": "day and night", "solar_day_hours": 24.0}, "planet")["bands"] == []


# --- the page ----------------------------------------------------------------------------------------------------------------------

@pytest.fixture
def with_row(client, monkeypatch):  # noqa: F811
    from planetgen.web.lib import apiclient

    def install(row):
        detail = pages._detail()
        detail["planets"][0].update(row)
        monkeypatch.setattr(apiclient, "get_system", lambda db, system_id: detail)
        return client
    return install


def test_a_planet_page_has_the_clock_the_note_and_the_bands(with_row):
    page = with_row(_row()).get(f"/object/{PLANET}").get_data(as_text=True)
    assert 'id="item-day-night"' in page and "data-day-counter" in page
    assert f'data-epoch-ms="{day_counter.EPOCH_MS}"' in page and 'data-day-hours="24.0"' in page and 'data-night-hours="12.0"' in page
    assert 'role="timer"' in page and 'aria-live="off"' in page, "a clock that moves must not chatter at a screen reader"
    assert "daycounter.js" in page
    assert "faster than the real rotation" in page
    assert '<th scope="row">Equator</th>' in page and 'class="day-counter-bands"' in page
    assert re.search(r"<caption>Day and night by latitude band</caption>", page)


def test_the_panel_sits_between_the_sensor_view_and_the_facts(with_row):
    page = with_row(_row()).get(f"/object/{PLANET}").get_data(as_text=True)
    assert page.index('id="item-day-night"') < page.index('id="item-facts"')


def test_a_locked_planet_page_says_so_and_has_no_clock_or_script(with_row):
    page = with_row(_row(day_night.LOCKED, **{f"band{n}_state": s for n, s in enumerate(RING_STATES, start=1)})).get(f"/object/{PLANET}").get_data(as_text=True)
    assert "Tidally locked" in page and "data-day-counter" not in page and "daycounter.js" not in page
    assert "Day centre" in page


def test_a_page_with_no_stored_state_has_no_panel(client):  # noqa: F811
    for ref in (PLANET, MOON, STAR):
        page = client.get(f"/object/{ref}").get_data(as_text=True)
        assert 'id="item-day-night"' not in page and "daycounter.js" not in page
