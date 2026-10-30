**Headline:** A system page now lists the orbits of its stars, so a triple or larger system shows what each orbit joins.

### Added
- **"Orbits of the stars" on the system page (MAP.213, slice C).** A system with two or more stars gets a table below its tree of bodies with one row per orbit in the stored tree of stars: what the orbit joins (two stars, or a star and the pair hung below it, named by its orbit number), whether it is a close, wide or compact pair, its mass, separation, period, eccentricity and how close the pair is (detached, semi-detached, contact and so on). A single star has no panel. `GET /api/systems/<id>` carries the same tree as `orbits` and each star's `parent_node_id` and `slot`; the API version does not change, since nothing was removed.
