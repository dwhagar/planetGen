**Headline:** The galaxy pages name the regions where sectors are charted instead of four quadrants.

### Changed
- **Quadrants are gone from the Galaxy page, the Sectors list and the sector page.** The Galaxy page's Quadrants table is now a Charted regions table (name, the family it belongs to, sectors charted and share charted, most charted first). The Sectors list shows each sector's region in place of its quadrant and no longer has a quadrant menu or a position sort. A sector page says which region it lies in. `GET /api/sectors` no longer takes `quadrant` or `facets` and each item carries its `region`. Systems still show their octant within the sector, which is a different thing.
