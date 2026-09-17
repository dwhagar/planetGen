# tests/test_wikiclient_mediawiki_integration.py

"""
End-to-end test for `wikiClient.mediawiki.MediaWikiBackend` against a real
MediaWiki instance -- confirms the login/token/edit request shape
`mediawiki.py` sends is one a real server actually accepts (the mocked
tests in `test_wikiclient_mediawiki.py` can only prove the backend behaves
correctly for a *given* response, not that MediaWiki would ever actually
send that response back).

Opt-in only, via the `mediawiki_config` fixture (`conftest.py`): skipped,
not failed, unless `PLANETGEN_TEST_MEDIAWIKI_BASE_URL`/
`PLANETGEN_TEST_MEDIAWIKI_USERNAME`/`PLANETGEN_TEST_MEDIAWIKI_PASSWORD`
point at a real, reachable MediaWiki instance whose bot password has
page-create rights. A disposable instance (e.g. MediaWiki's own official
Docker image) works well for this -- there is no cleanup step here to
delete the created page afterward, same reasoning as
`test_wikiclient_wikijs_integration.py`'s own note.
"""

import uuid

import pytest

from wikiClient import WikiClientPageExistsError
from wikiClient.mediawiki import MediaWikiBackend


def test_create_page_against_real_instance(mediawiki_config):
    base_url, username, password = mediawiki_config
    client = MediaWikiBackend(base_url, username, password)

    # Unique per run so repeated test runs against a long-lived instance
    # never collide with a page an earlier run left behind.
    title = f"PlanetgenTest/{uuid.uuid4().hex}"

    page = client.create_page(
        path=title,
        title=title,
        content="Created by planetGen's wikiClient mediawiki backend integration test.",
        summary="planetgen-test",
    )

    assert page.path == title
    assert page.id is not None
    assert page.url.startswith(base_url)


def test_create_page_duplicate_title_against_real_instance(mediawiki_config):
    base_url, username, password = mediawiki_config
    client = MediaWikiBackend(base_url, username, password)
    title = f"PlanetgenTest/{uuid.uuid4().hex}"

    client.create_page(path=title, title=title, content="First version.")
    with pytest.raises(WikiClientPageExistsError):
        client.create_page(path=title, title=title, content="Second version.")
