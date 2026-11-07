**Headline:** The log on the admin Update page and the log of a queue run now open scrolled to the last line and stay there as the log grows.

### Changed

- `src/html/static/logfollow.js` (UX.284): a log box marked `data-follow-end` is scrolled to its end when the page loads and whenever its text changes. A reader who scrolls up pauses the following, so the lines they are reading do not move; scrolling back to the end resumes it. The Update page reloads itself every few seconds while an update runs, so the pause and the place being read are kept in the browser's session storage across the reload (blocked storage means it follows the end). `/admin/update` and `/admin/queue/<run>` use it; the Generate job output already followed its end. Without script the box shows the log from its top, as before.
- `src/tests/test_logfollow_browser.py` tests the script in Chromium on a bare page.

generation-output: unchanged
