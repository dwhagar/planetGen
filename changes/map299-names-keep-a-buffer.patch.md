**Headline:** In a dense area no map name overlaps or crowds another, and where names must be cut the brighter or heavier thing keeps its name.

### Fixed
- **Names stood closer than their buffer, and sometimes over each other (MAP.299).** Three causes in the label engine: the screen grid only compared names in the same cell, so two names a few pixels apart across a grid line were never compared; a name's width was measured in one font whatever its kind, though a sector's name is semibold and a nebula's italic; and a name was placed at a fractional position and drawn at a rounded one. Names are now whole pixels wide enough for their own font, the grid reaches each name's buffer, and a browser test of the fixture's densest region finds none closer than the buffer at three zooms.
- **The buffer is a setting.** `map.labelBuffer` (pixels, 0 to 40, default 4) in the browser's preferences, read at each placement.
- **Where names of one kind compete for room the brighter or heavier keeps the better place.** A star by its brightness, an object or body by its mass, a nebula by its volume; the budget still chooses which things get a name and the order between kinds is unchanged.
