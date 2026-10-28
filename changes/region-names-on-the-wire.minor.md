**Headline:** The Galaxy Map's data now carries the name of every charted region, ready for the map to label.

### Added
- **Region names on the wire (API.35).** The node of a charted region (a 9-block holding at least one generated sector) has `name` ("Fardoo Beez Region") and `parent` (its family, "Fardoo Reach"), and so does the stage that lists such blocks as children and the stage of the region itself. Names come from the gated phoneme codec under the galaxy's naming key, so a new key renames them; uncharted blocks, other block sizes and a galaxy with no key carry none. The format document and API reference say so, and `tests/test_region_wire_names.py` pins the fields.

### Changed
- **The stored Galaxy Map nodes are rebuilt once** (`TILE_FORMAT` 5) so they pick up the names.
