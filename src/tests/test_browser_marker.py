# tests/test_browser_marker.py

"""
TEST.304: every test file that drives headless Chromium carries the `browser` marker, whatever it is called, so the
all-core phase of `scripts/test.py` never starts Chromium. The marker follows the file's own Playwright import or a
borrowed browser fixture (`conftest._uses_playwright`), so no file can import Playwright and miss it.
"""

import os
import types

from tests import conftest

TESTS = os.path.dirname(os.path.abspath(__file__))


def test_every_file_that_imports_playwright_is_a_browser_file():
    found = sorted(name for name in os.listdir(TESTS) if name.startswith("test_") and name.endswith(".py")
                   and conftest._imports_playwright(os.path.join(TESTS, name)))
    # The files the old name rule (test_web_a11y.py, test_web_browser_*) missed.
    assert {"test_architecture_map.py", "test_web_admin_menu_fits.py", "test_web_adminshell_browser.py",
            "test_web_edit_menu.py", "test_web_generate_alignment.py", "test_web_generate_flow_browser.py",
            "test_web_header_menus.py", "test_web_job_terminal.py", "test_web_settings_browser.py",
            "test_web_system_nav_row.py", "test_web_ui_fields_browser.py", "test_web_a11y.py"} <= set(found)
    assert "test_browser_marker.py" not in found and "test_fuzz_web_routes.py" not in found


def test_the_import_forms_and_a_borrowed_fixture_count():
    name = "play" + "wright"  # spelt apart, so this file's own text imports nothing
    for text in (f'pytest.importorskip("{name}.sync_api")', f"sync_api = pytest.importorskip('{name}.sync_api')",
                 f"from {name}.sync_api import sync_{name}", f"import {name}"):
        assert conftest.PLAYWRIGHT_IMPORT.search(text), text
    assert not conftest.PLAYWRIGHT_IMPORT.search(f"# the {name} fixtures live in test_web_a11y")
    here = os.path.abspath(__file__)
    assert conftest._uses_playwright(types.SimpleNamespace(fixturenames=("gl_browser", "page"), path=here))
    assert not conftest._uses_playwright(types.SimpleNamespace(fixturenames=("page",), path=here))
