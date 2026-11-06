**Headline:** Every body an asteroid belt keeps one by one now has a page of its own, with its mass, size, makeup, orbit and the history of its belt, and the belt's list of bodies links to it.

### Added
- **A page for each stored body of a belt (GEN.1070).** `/object/belt:<ID>/member/<rank>` shows the body by its designation (it has no name), its shape word, mass, diameter, density, type and makeup, the orbit it was saved with (mean, closest and farthest distance, year, eccentricity, inclination and the three angles), its share of the belt's mass, and where its belt came from and what shaped it, with links to the heavier and the lighter neighbor, the belt and the belt's list. The list on the belt's information page links each designation to its page.
- **`GET /api/objects/<ref>/belt-members/<rank>`** answers the same data. It is a new route, so the API version does not change. A rank above the number of bodies the belt keeps is drawn on demand and has no page yet; it answers 404, as does a rank not written the way the designation writes it.

generation-output: unchanged

TODO items resolved: none (GEN.1070 stays open for the designation display and the dots on the system map).
