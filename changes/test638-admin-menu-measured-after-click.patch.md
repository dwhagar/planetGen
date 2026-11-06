**Headline:** The admin menu layout test measures the button and its panel at the same moment, after the click, so it no longer fails on pages whose map scrolls the page on load.

### Fixed
- `test_web_admin_menu_fits` failed for the laptop, phone, short and tablet windows on the sector page: the admin panel was "not at its button", the button at a negative y. The page was not wrong. A page whose map fills its window scrolls itself to the map when it has loaded (MAP.350, `mapfit.js`; the sector page scrolled 291 px at the laptop size), so the button lay 73 to 235 px above the window when the test measured it, before the click that brings it back into view, and the panel was measured after. The test now measures the Admin button, and the Navigate button, after the click, with the panel. The check itself is unchanged: the panel must hang just below the button or stand just above it, and inside the window.

generation-output: unchanged
TODO items resolved: TEST.638
