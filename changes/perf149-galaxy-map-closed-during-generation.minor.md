**Headline:** The Galaxy Map is closed while a galaxy generation runs.

### Changed
- **No map from a half-built galaxy.** While a run is changing the galaxy, the Galaxy Map page says a generation is running instead of drawing the map: it shows the stage the run is at and the banner's sentence, checks again by itself every 30 seconds, and shows the map as soon as the run has ended. An admin gets a link to the Queue page. The page answers 503 and is not kept by browsers or caches.
