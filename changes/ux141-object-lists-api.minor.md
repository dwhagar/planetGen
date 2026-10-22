**Headline:** Every kind of object has a list that reads the object index, a page at a time without skipping rows.

### Changed
- **`GET /api/object-lists/<kind>` (UX.141, first part).** One page of the systems, stars, planets, moons, belts, comets, facilities or any phenomenon type from the object index, sorted by name, class, mass or radius, filtered by a name prefix, a class or a sector. The next page starts after the last row's object ID, so a deep page costs the same as the first and no query skips rows with `OFFSET`. The index is filled by the warm-up (`planetgen warm object-index`).

TODO items resolved / Added: none / Added: none
