### Fixed

- Galaxy Map: moving the view (right-drag or Shift-drag) now stops one and a half views from where the stage opened, as intended; before, nothing held it and the view could slide away for good.
- Galaxy Map: reloading the page after the map's Back button keeps Forward working.

### Added

- Tests for the page scripts under node (`src/tests/js`, run by `test_js_unit.py`): the Galaxy Map's drill-down, history, address bar, zoom, pan and tilt limits and buttons; the Sector Map's zoom, turning and picking; the phenomenon diagram's zoom; the Generate page's job panel; the facility form (TEST.57, TEST.58).
- Browser tests: every map button changes the view, the System Map's selection and measuring, the Galaxy Map drill-down by clicks with Back and Forward and the free camera (TEST.55, TEST.58, TEST.59), and no overlapping or off-screen controls on any page at 390, 600, 820 and 1280 px in both themes (TEST.56). The large-nebula diagram's dead "-" button is pinned as a known failure for UX.21.
