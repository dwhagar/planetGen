**Headline:** The star, planet and moon lists filter and sort by the same stored physical values as the search page: gravity, temperature, pressure, orbit, dose, luminosity, age, gases and the rest.

### Added
- `/lists/star`, `/lists/planet` and `/lists/moon` carry a "Filter and order by physical values" form with the search page's range boxes and unit menus (`planet_min_gravity=0.8&planet_max_gravity=1.2`, `planet_value_sort=-gravity`), chips that take one range away, and a Values column for the values a range or an order is on that no column shows already. The columns of a physical value (orbit, surface temperature, a star's temperature, luminosity and age) sort the list by it. The table keeps the ranges while it is sorted and paged, and a submitted form is redirected to its short address like the search page's. The partial `value_group.html` is shared with the search page.
- `GET /api/object-lists/<kind>` for those three kinds takes `sort=<key>` for any numeric key of the search and `<kind>_min_<key>`, `<kind>_max_<key>`, `<kind>_<key>_unit` range filters; items carry `values` for the keys in use. Value sorts and ranges read the kind's own table through a join on the object index, keyset paged.

### Changed
- A sort the star list did not know before (`age`) is now a valid sort rather than a 400.

generation-output: unchanged
TODO items resolved: UX.228 part 1 (the star, planet and moon lists); the Systems list and Nearby are the other parts of the item
