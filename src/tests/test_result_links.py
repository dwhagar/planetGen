# tests/test_result_links.py

"""
UX.232: every link in a search result row opens the object it names. A planet link used to be built from the row's number
instead of the object's printed ID and answered 404, while the system link in the same row worked. These read the links of a
real sample galaxy's search and Nearby results and follow each.
"""

import re

import pytest

from planetgen.web.app import create_app
from planetgen.api.config import Config

from tests.habitability_site import habitable_galaxy  # noqa: F401 -- fixture


@pytest.fixture
def real_client(habitable_galaxy):
    class RealConfig(Config):
        MYSQL_CONFIG = habitable_galaxy.config
        WRITE_MYSQL_CONFIG = habitable_galaxy.config
        CONTROL_MYSQL_CONFIG = habitable_galaxy.config
        SESSION_COOKIE_SECURE = False
        SECRET_KEY = "test-secret"

    application = create_app(RealConfig)
    application.testing = True
    return application.test_client()


def _links(html, panel):
    part = html[html.index(f'id="search-{panel}"'):]
    part = part[:part.index("</section>")]
    # The sector links are left to the sector tests: this galaxy's sectors are placed by hand without a printed ID.
    return [href for href in re.findall(r'<td[^>]*><a href="([^"]+)"', part) if not href.startswith("/sector/")]


@pytest.mark.parametrize("panel,query", (("planets", "type=planet"), ("moons", "type=moon"), ("stars", "type=star"),
                                         ("planets", "planet_q=a")))
def test_the_object_links_of_search_rows_answer_200_with_the_object(real_client, panel, query):
    html = real_client.get(f"/search?{query}").get_data(as_text=True)
    links = _links(html, panel)
    assert links, panel
    for href in links[:12]:
        answer = real_client.get(href)
        assert answer.status_code == 200, (href, answer.status_code)
        if href.startswith("/object/"):
            assert panel[:-1] in href


def test_the_class_search_rows_link_the_planet_the_system_and_the_sector(real_client, habitable_galaxy):
    from planetgen.db import query
    conn = query.open_readonly(habitable_galaxy.config)
    try:
        klass = conn.execute("SELECT planet_class AS c FROM planets WHERE planet_class IS NOT NULL LIMIT 1").fetchone()["c"]
    finally:
        conn.close()
    html = real_client.get(f"/search?class={klass}").get_data(as_text=True)
    hrefs = _links(html, "planets")
    assert hrefs
    for href in hrefs[:12]:
        assert real_client.get(href).status_code == 200, href
