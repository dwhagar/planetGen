# tests/test_web_formulas.py

"""
UX.149 on the site: the `formula` and `chem` filters, and no page of the sample galaxy (every page the accessibility
test loads, and the class, species, phenomenon and system pages they link to) shows a formula that is not written
out: what is left of a page's text once its tags are gone holds no formula `planetgen.util.chem.split` would mark.
Skipped without the MySQL test server, like `test_web_a11y.py`.
"""

import html
import re

from planetgen.api.authz import SESSION_COOKIE_NAME
from planetgen.util import chem

from tests.test_web_a11y import (  # noqa: F401 -- fixtures
    admin_token,
    page_targets,
    sample_job,
    sample_job_tree,
    sample_params,
    site_app,
    site_db,
)

LINKED = re.compile(r'href="(/(?:classes|species|phenomenon|system)/[^"#?]*)"')


def _visible_text(page):
    """A page's text as a reader sees it, each tag a space (so `H<sub>2</sub>O` reads `H 2 O`, no formula)."""
    page = re.sub(r"(?s)<(script|style)\b.*?</\1>", " ", page)
    return html.unescape(re.sub(r"(?s)<[^>]+>", " ", page))


def test_the_filters_write_a_formula_out(site_app):
    with site_app.test_request_context("/"):
        page = site_app.jinja_env.from_string("{{ 'H2O'|formula }} {{ 'H+, e- and <b>'|chem }} {{ ''|chem }}").render()
    assert page == "H<sub>2</sub>O H<sup>+</sup>, e<sup>−</sup> and &lt;b&gt; "


def test_no_page_shows_a_formula_that_is_not_written_out(site_app, page_targets, admin_token):
    client = site_app.test_client()
    client.set_cookie(SESSION_COOKIE_NAME, admin_token)
    paths = [path for path, _admin in page_targets.values() if path]
    linked = set()
    for path in paths:
        linked |= set(LINKED.findall(client.get(path).get_data(as_text=True)))
    found = []
    for path in paths + sorted(linked - set(paths)):
        response = client.get(path)
        if response.status_code != 200 or response.mimetype != "text/html":
            continue
        page = response.get_data(as_text=True)
        if path == "/classes/nebula/M":
            assert "Dark nebula: H<sub>2</sub>, He, CO, PAHs, silicates." in page
        found += [(path, part) for part in chem.split(_visible_text(page)) if isinstance(part, chem.Formula)]
    assert "/classes/nebula/M" in linked
    assert not found, found
