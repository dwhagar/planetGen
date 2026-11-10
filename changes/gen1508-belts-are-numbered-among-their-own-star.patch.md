**Headline:** A belt is now numbered among the belts of its own star, closest first from 1, so a lone belt behind four planets is Asteroid belt 1, not Asteroid belt 5.

### Changed
- **Belt numbers count belts (GEN.1508, Boss's live observation).** A belt has no stored name; the pages, the object lists, the nearby list, the parent chain and the object index built "Asteroid belt N" from the belt's `orbital_index + 1`, which is its slot among all the star's bodies (planets and belts together). They now count the belts of the star by their place out from it, closest first, from 1 (`names.bodies.belt_number`; the database side is `db/belt_names.py`, and the object index's belt name is `objectindex.BELT_NAME_SQL`). Each star of a wide pair counts its own belts, only belts that exist count (no gaps), and a lone belt is `Asteroid belt 1`. Comets carry stored names and no number, and no dwarf planet is named by slot, so they have no such fault. Nothing generated or stored in the body tables changes; the `belt` rows of the derived `object_index` table get the new names the next time it is built (`objectindex build`, which rewrites only changed rows).

generation-output: unchanged
generation-affects: none
