# tests/test_wikiclient_wikijs_integration.py

"""
End-to-end test for `wikiClient.wikijs.WikiJsBackend` against a real
Wiki.js instance -- confirms the GraphQL request shape `wikijs.py` sends is
one a real server actually accepts (the mocked tests in
`test_wikiclient_wikijs.py` can only prove the backend behaves correctly
for a *given* response, not that Wiki.js would ever actually send that
response back).

Opt-in only, via the `wikijs_config` fixture (`conftest.py`): skipped, not
failed, unless `PLANETGEN_TEST_WIKIJS_BASE_URL`/`PLANETGEN_TEST_WIKIJS_TOKEN`
point at a real, reachable Wiki.js instance whose token has page-create
rights. A disposable instance (e.g. Wiki.js's own official Docker image)
works well for this -- there is no cleanup step here to delete the created
page afterward (Wiki.js's GraphQL API supports `pages.delete`, but adding
that isn't worth it just for a test fixture against a throwaway instance
meant to be discarded after the run).
"""

import uuid

import pytest

from wikiClient import WikiClientPageExistsError
from wikiClient.wikijs import WikiJsBackend


def test_create_page_against_real_instance(wikijs_config):
    base_url, api_token = wikijs_config
    client = WikiJsBackend(base_url, api_token)

    # Unique per run so repeated test runs against a long-lived instance
    # never collide with a page an earlier run left behind.
    path = f"planetgen-test/{uuid.uuid4().hex}"

    page = client.create_page(
        path=path,
        title="planetGen wikiClient wikijs backend test",
        content="# Test page\n\nCreated by planetGen's wikiClient wikijs backend integration test.",
        description="Throwaway page created by an automated test.",
        tags=["planetgen-test"],
    )

    assert page.path.strip("/") == path
    assert page.id is not None
    assert page.url.startswith(base_url)


def test_create_page_duplicate_path_against_real_instance(wikijs_config):
    base_url, api_token = wikijs_config
    client = WikiJsBackend(base_url, api_token)
    path = f"planetgen-test/{uuid.uuid4().hex}"

    client.create_page(path=path, title="First", content="First version.")
    with pytest.raises(WikiClientPageExistsError):
        client.create_page(path=path, title="Second", content="Second version.")
