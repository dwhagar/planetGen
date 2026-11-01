**Headline:** A planet or moon's page now has a Day and night panel: a clock that counts its own day at the real rate, the state of its day and night in words, and the day and night hours of each latitude band beside the band's temperature.

### Added
- **The day-and-night counter (UX.243).** On the page of a planet or moon that has a stored state of day and night, a clock counts the object's own solar day at its real rate: the day number since 1 January 2000, 12:00 UTC (the same moment for every object, so two pages agree), the time into the day, and whether it is day or night with the time left of it. A tidally locked body has no clock; the panel says that one face stays toward its star (or its planet, for a moon), with a permanent day side and night side. A body that has no night says so. The page says that the picture of the body turns much faster than the real rotation, so the picture's speed and the clock stay apart. Beside the clock the seven latitude bands (or the seven rings of a locked world) are laid out with their annual mean temperature, their day and night hours and their state, each number with its unit. The values are the ones already stored (GEN.1097, GEN.1108), so there is no schema or API change; a body made before them has no panel.

### Fixed
- **The browser test's bright-star rows.** `test_web_a11y.py` wrote bright-star rows one column short of what the scatter writes now, so every page check failed building its database; the row has its last column again.

generation-output: unchanged

TODO items resolved: UX.243.
