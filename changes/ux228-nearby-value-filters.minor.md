**Headline:** What's nearby filters and orders what it lists by the same stored physical values as the search page: planets within 20 parsecs with 0.8 to 1.2 g, stars by age, moons by orbit.

### Added
- `/nearby` carries the search page's "Filter and order by physical values" form (`planet_min_gravity=0.8&planet_max_gravity=1.2`, `star_value_sort=-age`): a range lists only the stars, planets or moons in the sphere whose stored value lies in it, a kind with a range is listed even when it is not ticked, an order by a value lists that kind's bodies by it (nearest first among equals; planets, then moons, then stars when several are set). Each range shows as a chip with its unit, the rows carry a Values column, and a submitted form is redirected to its short address.
- `GET /api/near` takes the same ranges (`<entity>_min_<key>`, `<entity>_max_<key>`, `<entity>_<key>_unit`) and `<entity>_sort=<key>` / `-<key>`; rows carry `values`. `nearest=K` takes the ranges; an order by a value with `nearest`, a range on another kind than the one ordered, and a key a kind does not have are a 400. Value orders hold one heap key per candidate body of the sphere and name only the page's rows.

generation-output: unchanged
TODO items resolved: UX.228 part 3 (Nearby); the Systems list is the remaining part of the item
