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


def test_every_browser_fixture_that_starts_playwright_stops_it():
    """TEST.696: `test_life_text_a11y_browser.py` started Playwright and never stopped it, so its asyncio loop stayed alive for every later module in the same xdist
    worker and about 1,500 browser tests failed to start. A file that calls `.start()` on `sync_playwright()` must call `.stop()` in the `finally` after its `yield`; one that uses it as a
    context manager is closed by that."""
    import pathlib
    import re

    folder = pathlib.Path(__file__).resolve().parent
    starters = []
    for path in sorted(folder.glob("test_*.py")):
        text = path.read_text(encoding="utf-8")
        if re.search(r"sync_playwright\(\)\s*\.start\(\)", text) and not re.search(r"yield[^\n]*\n\s*finally:[\s\S]{0,300}?\.stop\(\)", text):
            starters.append(path.name)
    assert not starters, f"these files start Playwright and never stop it: {starters}"
