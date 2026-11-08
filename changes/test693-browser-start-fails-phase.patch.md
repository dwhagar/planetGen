**Headline:** A browser that will not start now fails the test run instead of skipping 1,500 tests.

### Fixed
- **Browser tests no longer skip silently (TEST.693).** The test preflight demanded Playwright 1.56.0 whatever Chromium the machine carried, so on an image with Chromium 1243 every browser test skipped with "Executable doesn't exist" and the phase still printed a pass. The preflight now asks the installed Playwright for its own Chromium and, if it is not there, names the Playwright that matches each Chromium build on the machine; the cloud setup script picks the matching release too. In the browser phase of `scripts/test.py` a test that skips because Playwright or Chromium would not start is now a failure.

generation-output: unchanged

TODO items resolved: TEST.693.
