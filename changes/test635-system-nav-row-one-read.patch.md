**Headline:** The test that the system page's buttons share one row reads them in one step, so the page's own scroll to the System Map no longer makes it fail.

### Fixed
- TEST.635: `test_system_buttons_share_one_row` read each button's box with a call of its own, and the System Map scrolls the page to itself a moment after the page loads (MAP.350), so the boxes were read either side of that scroll (one at y 217, the next at y -170) and looked like two rows at every width. The three boxes are now read in one call. The Nav menu test also expects Add waypoint as the last item (NAV.49).
