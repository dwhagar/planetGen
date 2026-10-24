**Headline:** The species tech filter test no longer fails when the population pass makes a species above tech 6.

### Fixed
- **`test_api_species_tech_filters` failed about once in two full runs (TEST.299).** The test raised the first species to tech 6.0 and expected it to lead the tech-descending list, but the population pass draws each species' tech from the galaxy fixture, and now and then makes one higher (6.55 in the reported run), which then led the list. The test now sets every species to tech 1.0 first, so the one it raises is the top.
TODO items resolved: TEST.299. Added: none.
