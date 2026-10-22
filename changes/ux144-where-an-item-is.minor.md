**Headline:** Every item's page now says where it is, and a nebula lists what is inside it.

### Changed

- A star, planet, moon, belt or comet page, a system page and a phenomenon page each have a Location block: the places above the item as links (galaxy, sector, system, the planet a moon circles) and a line `Inside: <nebula or supernova remnant>` linking the cloud it sits in, or a line saying it sits in none.
- A nebula's or supernova remnant's page lists what is inside it, grouped by kind (star systems, rogue planets, black holes and the rest) with a link to each, or says nothing generated sits inside it.
- New API route `GET /api/phenomena/<type>/<id>/contents` behind that list.

TODO items resolved / Added: UX.144 (nebulae and supernova remnants; globular clusters and superbubbles are not stored as containers yet) / Added: none
