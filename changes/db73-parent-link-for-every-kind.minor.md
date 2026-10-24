**Headline:** Every body now records what it orbits in one way, so a planet can leave its star and keep its identity.

### Changed
- **One parent link for every kind (DB.73).** Stars, asteroid belts, comets, black holes, neutron stars, rogue planets and interstellar comets now carry the `parent_type` and `parent_id` that planets and moons already had (a star, a planet, a system, an orbit node, a rogue planet, or free), filled in for everything already stored. A moon can be stored under a rogue (free) planet, with its own ID under that planet's ID, and removing the rogue planet removes its moons. Changing a planet's parent rewrites only the link: the row and the ID stay. Database version 140 (migration 0140); run `./update.sh`.
