**Headline:** A moon of a free (rogue) planet can now be found by name, looked up by its ID and opened on its own page.

### Changed
- **Moons of free planets are read where other moons are (DB.94).** The object index, the object lookup (`GET /api/objects/<ref>`), the Search page's Moons panel and the moon's own page all handle a moon that orbits a free planet instead of a star: the lookup names the rogue planet as its parent and gives its galaxy point, Search lists it with the free planet where a planet's name goes (and no system), and its page links to the free planet in place of a system. Nothing generates such moons yet.
