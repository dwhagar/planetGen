**Headline:** The Uncharted stars and Phenomena tables page by key, so a page far down the list costs what the first page does.

### Changed
- **Uncharted stars page by key (PERF.75).** `GET /api/uncharted-systems` takes `after=<star id>` and the table passes the last star of each page on, so the brightest-first list and the sector sort walk their index to the page instead of skipping rows. Equal values are now ordered by id in the sort's own direction, which lets a descending page read the index backward with no sort.
- **The Phenomena table pages the scattered tail by key.** `GET /api/phenomena` takes `after=<scatter id>`; a page that ends among the unbuilt scattered phenomena hands the next page its last id.
- **The Systems name sort reads either way round from its index.** Systems with the same name are ordered by id in the sort's direction, so the descending name sort needs no sort of the table either.
