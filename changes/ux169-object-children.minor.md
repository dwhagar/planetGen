**Headline:** `GET /api/objects/<ref>/children` now pages through what belongs to a sector, system, star or planet, with the counts a tree row needs.

### Added

- **The children of an object (UX.169).** `GET /api/objects/<ref>/children` answers one keyset page of what is inside a sector (its systems, then the phenomena placed in it), a system (its stars, then the planets, belts and comets that orbit the whole system), a star (its planets and belts in their order out from it, then its comets) or a planet or free planet (its moons). Each item carries a short summary of its own row, the number of stars, planets, moons, belts and comets anywhere below it and how many it holds directly, so a tree row can say "3 planets, 5 moons" and show its expand arrow without another request. A sector with thousands of systems is read 200 at a time by the previous page's token, never by offset. This is the feed the summary trees of the coming object pages use.
