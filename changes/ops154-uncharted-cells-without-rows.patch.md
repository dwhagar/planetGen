**Headline:** The Sectors list, its counts and the sector page now find an uncharted cell that holds a waiting star or phenomenon from its address and the block light parts, so such a cell no longer needs a sector row of its own (part 1 of OPS.154; the rows stay until part 2 stops making and deletes them).

### Changed
- **Cells without a row are listed (OPS.154).** With `coverage=uncharted|all`, `list_sectors` goes on, after the rows, into the cells that hold a part of the block lights, have no chart-ledger entry and no sector row (`query.virtual_cells`): by distance (shell, band, wedge; backward for a descending sort), read by key so any page costs the same, with `after` and `offset` both working across the rows and the cells. Such an item has `virtual: true` and its sector ID (the designation) as `id`. The `quadrant` filter works from the cell's position; the `systems` and `density` sorts count a cell as zero. A cell that has a row (a rogue's) is listed once, as its row.
- **Counts.** `count_uncharted_sectors` and `count_sectors` for `uncharted` and `all` add `count_virtual_cells`; the placeholder-row count is kept by `countcache` now, as it is read on every list.
- **The cell answers by its ID.** `GET /api/sectors/<id>` for a designation no row has answers for a cell that holds a waiting object (`query.virtual_sector_detail`), so the web sector page opens it as the uncharted page; `after=<id>` takes such an ID.
- **While the rows are still there** (more than `VIRTUAL_PLACEHOLDER_LIMIT`, 100,000, of them) the cells are not listed from the parts, since they have rows and the walk would read the whole table: the list is the rows alone, as before.
- Not in this part: search by name does not find a cell with no row (a designation typed in the box) until part 2.

generation-output: unchanged
generation-affects: none
