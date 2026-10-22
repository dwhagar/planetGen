**Headline:** The cache warm-up also builds the Lists page and the first page of every kind's list.

### Changed
- **Warm the new list pages (PERF.141).** The `lists` warm family now visits `/lists` and `/lists/<kind>` for every kind, so the first visit to a list of stars, planets, moons, belts, comets or a phenomenon type is answered from the cache. Deeper pages are not warmed, because a page of these lists starts after the last object of the one before and costs the same at any depth.

TODO items resolved / Added: none / Added: none
