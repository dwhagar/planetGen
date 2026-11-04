**Headline:** The nav chart and System Map box-fit script loads its sibling the way every other page script does, so the page-shell test passes again.

### Changed
- **TEST.568.** `static/boxfit.js` imported `mapfit.js` with a static `import`, which `test_page_shell` refuses because a static import drops the page's version query and could load a stale sibling. It now reads the version query from its own URL and imports `mapfit.js` dynamically with it, as the other page scripts do. No generator output changes.

generation-output: unchanged

TODO items resolved: TEST.568.
