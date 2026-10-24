**Headline:** Every test that drives headless Chromium now runs in the browser phase of the test runner, not beside the all-core tests.

### Changed
- **Eleven Playwright test files were not marked `browser` (TEST.304).** `conftest.py` marked only `test_web_a11y.py` and `test_web_browser_*.py`, so `test_architecture_map`, `test_web_admin_menu_fits`, `test_web_adminshell_browser`, `test_web_edit_menu`, `test_web_generate_alignment`, `test_web_generate_flow_browser`, `test_web_header_menus`, `test_web_job_terminal`, `test_web_settings_browser`, `test_web_system_nav_row` and `test_web_ui_fields_browser` started Chromium in the all-core phase. A test is now marked by what it does: its file imports Playwright, or it takes a browser fixture (`conftest._uses_playwright`), so a new browser file cannot miss the marker. Test in `test_browser_marker.py`.
TODO items resolved: TEST.304. Added: none.
