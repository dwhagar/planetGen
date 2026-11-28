**Headline:** Phenomena: the core of a remnant saved with no sector has an ID of its own

### Fixed

- NAV.81: a supernova remnant saved with no sector (`planetgen generate supernova-remnant` with no sector) got an Object ID for the remnant but none for its collapsed core, so `/phenomena` listed the core by its row number and its page answered 404 for that number. The core now takes its ID with the remnant (`assign_uids`), keeping the name `<remnant> Core`, so the list and the page agree.
