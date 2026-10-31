**Headline:** Superbubbles are off until the Menu's button shows them, the Key is a small corner box, every charted region is a point on the galaxy view with a link that zooms to it, the sector card's buttons stay in view, the sector page and 3D system view fill the window like the Galaxy Map, and an opened sector stands alone, close and full of systems.

### Fixed
- **The Key covered too much of the map (MAP.356).** The Key is one column of short lines in a box no wider than two fifths of the map (13 rem at most), in the map's lower left corner.
- **Charted sectors and regions could not be found from the galaxy (MAP.357).** From the galactic view every filled region far enough from the others on the screen to be told apart gets its own point of interest with a zoom-in link that frames the whole region; regions close together share one point. A blank Go to box lists the five regions with the most generated sectors (`Charted places`), each a link.
- **The sector card's View system and Expand buttons could be out of view (MAP.358).** The card keeps at least 13 rem of height inside the Galaxy Map, so its buttons stay in the card and the map on a small window.
- **The sector page's map and the 3D system view were smaller than the Galaxy Map (MAP.359).** Both now fill 95 percent of the browser window height less a held header, the same fit the Galaxy Map uses.
- **An opened sector shared the picture with every other sector (MAP.360).** While a sector is opened the other sectors' blocks fade to a faint level so only the sector being viewed stands out, and the camera fits it with a 2 percent margin instead of 10.

### Changed
- **Superbubbles are off by default (MAP.362).** The map opens without them; the Menu has a Superbubbles button (and the address `?show=superbubbles`) that draws them, and the Key names them only while they are shown.
- **An opened sector shows more systems (MAP.361).** The marker budget in an opened sector is 1500 instead of 240 and marks merge only within 12 px instead of 22 px, so the view holds every system that fits and thins only when the window cannot.

No generator output changes.

generation-output: unchanged

TODO items resolved: none (TODO retires MAP.356 to MAP.362 on merge).
