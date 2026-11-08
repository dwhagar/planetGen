**Headline:** One browser test file that never stopped Playwright no longer breaks the 1,500 browser tests after it.

### Fixed
- **The browser phase no longer collapses after one file (TEST.696).** `test_life_text_a11y_browser.py` started Playwright for its tests and never stopped it. Playwright's event loop then stayed alive in that test worker, and the next browser file to start Playwright there failed with "Sync API inside the asyncio loop", which in a full run took about 1,550 browser tests with it. The fixture now stops Playwright, and a test checks that every file that starts Playwright stops it.

generation-output: unchanged

TODO items resolved: TEST.696 (the loop; the remaining browser failures follow in their own notes).
