**Headline:** Planets and moons on the System Map now show their ocean, land and ice, drawn from the stored water fractions and temperature.

### Added
- **The body look record (MAP.271).** Every star, planet and moon marker on the System Map carries one `data-look` JSON record built by `planetgen/web/maps/bodylook.py` from the stored row: the body family, a seed from its id, a base colour, the water fractions, ocean class and colour, hydrosphere, pressure, temperature and atmosphere rim colour. A value that is not stored is left out, not zero. It replaces the loose `data-color`, `data-hasatmosphere` and the temperature text the page used to parse.
- **The shared painter (MAP.272).** `static/bodytexture.js` holds the seeded noise, the 128-body texture cache, a paint queue that paints one texture per frame (largest first) and the size tiers (point, disc, 128 x 64, 256 x 128).
- **The solid-world painter (MAP.273).** `static/solidpaint.js` paints a planet or moon as a rough few-tone picture: the sea covers the stored ocean fraction of the sphere, ice sits where the temperature is below freezing (polar caps on a temperate world, all ice on a frozen one), airless worlds show craters. The info panel's Colour line says it is an estimate from temperature and stored water fractions until the composition colour (MAP.274) replaces it. Gas giants and stars keep their current look for now (MAP.275, MAP.276).

### Changed
- **A planet's sphere is no longer coloured by its class letter.** The flat marker keeps the class colour as the key to the class letter; the sphere takes its colour from the look record.
