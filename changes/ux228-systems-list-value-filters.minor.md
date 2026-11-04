**Headline:** The Systems list filters by what a system holds: systems with a planet at 0.8 to 1.2 g, a star of a given age, or a moon in a given orbit.

### Added
- `/systems` carries a "Filter by what the system holds" form with the search page's range boxes and unit menus for stars, planets, and moons (`planet_min_gravity=0.8&planet_max_gravity=1.2`). A range keeps the systems that hold at least one body of that kind whose stored value lies in it; ranges on more than one kind all have to be met. Ranges show as chips that take one away, the table keeps them while it is sorted and paged, and a submitted form is redirected to its short address.
- `GET /api/systems` takes the same ranges (`<entity>_min_<key>`, `<entity>_max_<key>`, `<entity>_<key>_unit`); `total` and `facets` follow them. `db.query.list_systems`, `count_systems` and `systems_facets` take `values`.

### Changed
- The value chips of the search page, the lists, and What's nearby come from one function (`searchpage.value_chips`); their text is unchanged.

generation-output: unchanged
TODO items resolved: UX.228 (part 2 of 3: the Systems list; the star, planet, and moon lists and Nearby merged earlier)
