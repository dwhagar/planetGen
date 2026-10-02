### Fixed
- **Faint stars are brighter on the Sector Map and the Galaxy Map (MAP.87).** One shared curve draws the dimmest red dwarfs with four times the halo light, tapering to no change at 1000 L☉ and up, slowly enough that a brighter star is never drawn fainter than a dimmer one (a Sun gets about 2.2 times). Display only: nothing stored changes.
- **Rogue planets on the Sector Map are faint until marked (MAP.82, MAP.84).** Unmarked, a rogue planet is a dim speck with no glow, smaller than any star and only picked by a click right on it; "Mark rogue planets" makes each one bigger, fully lit, glowing and ringed, with a wide pick reach.
- **"Mark rogue planets" shows when it is on (MAP.83).** It starts off and stays highlighted while on, in both themes.
- **No link out of a course being picked (NAV.30).** While choosing a NAV start or destination, the Sector Map's info panel offers only the pick button, and the Galaxy Map's star and cloud panels drop "View system" and "View phenomenon".
- **Galaxy Map bookmark keys no longer clash with the browser's tab keys (MAP.81).** The first nine bookmarks open with plain 1 to 9 while the map has focus, instead of Ctrl+1 to Ctrl+9.

### Added
- **Browser tests for the maps that need no database (TEST.70).** `test_web_browser_fixture_maps.py` drives the Galaxy Map and the Sector Map, served by the real Flask views over fixture data, through picking, hover, keys, Back and Forward, URL state, bookmarks and the scale line.
- **One module for the maps' shared helpers (MAP.63).** `static/mapcore.js` holds what the Galaxy Map, the Sector Map and the System Map each had their own copy of: reading the scene data, theme colors, info-panel fields, the highlight ring, the scale bar's numbers, fitting the canvas and picking a point of light on screen. No visible change.
- **One camera and input controller for the Galaxy Map and the Sector Map (MAP.64).** `static/mapcontrol.js` turns, moves and zooms both maps, with a zoom policy each view sets (free, a short range, or locked), and tells a drag from a click in one place. No visible change.
