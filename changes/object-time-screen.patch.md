**Headline:** The admin Stats page now shows how long one object of each kind took to generate, which the generator has been recording without showing.

### Changed

- The Generation speed panel of the Stats page gains a "Seconds per object" table: for each kind of generated object (system, star, planet, moon, belt ...) the objects timed, the mean, fastest and slowest time and its share of all the object time, over every saved run, longest-running kinds first. `GET /api/admin/generation-stats` returns it as `object_types`, now a list of rows instead of a kind-to-mean map. `object_timing.kind_summaries` replaces `type_means` and the unused `run_totals`. generation-output: unchanged
