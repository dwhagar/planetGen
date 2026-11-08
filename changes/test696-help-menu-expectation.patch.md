**Headline:** The header menu test reads the Help menu's list from the help registry instead of a copy that went stale.

### Fixed
- **Help menu test (TEST.696).** `test_web_header_menus.py` expected the Help menu to hold a fixed list that no longer matched it (Properties had left, The sky behind a system had joined), so ten of its cases failed. It now reads the expected list from the same registry the menu is built from.

generation-output: unchanged

TODO items resolved: TEST.696 (header menus).
