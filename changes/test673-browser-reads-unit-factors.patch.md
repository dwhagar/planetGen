**Headline:** The browser's sky and day counter read the parsec, the year and the hour from the unit registry instead of writing the numbers themselves.

### Changed

- **Unit factors read from the registry in the browser (TEST.673).** `skysphere.js` took the parsec in kilometers and the year in seconds as numbers written in the file, and `daycounter.js` wrote 3600 in its clock text. Both now get them from `unitregistry.js` (`convert`), which is what `test_no_browser_script_declares_a_unit_factor` requires; the values are the same, so the sky and the clock do not change.

generation-output: unchanged
