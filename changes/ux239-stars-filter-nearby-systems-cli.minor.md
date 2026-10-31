**Headline:** The stars filter now also works in Nearby, the Systems list and its API, and the command line, and a system can be sorted by how many stars it holds.

### Added
- **The stars filter everywhere a system is listed (UX.239, the later slices).** What's nearby has a Stars menu (Single, Binary, Trinary, Quaternary, Quintenary and beyond) that keeps the systems with that many stars and what is in them, and a Stars column on its rows; `GET /api/near` takes `stars=` (a count, a count and a plus, or a name, as in `/api/search`) and a `system` row carries `star_count`. `GET /api/systems` takes `stars=` too and a `sort=stars` that orders by the number of stars, and the Stars column of the Systems tables sorts by it. The command line takes `--stars` on `query systems` and `query near` (the help names every choice) and says how many stars each system holds (single, binary, trinary, quaternary, or the count). A value that is none of these is a 400 from the API and an argument error on the command line; a hand-edited Nearby address with a bad value is ignored. Phenomena and facilities outside a system are left out of a Nearby list with the filter on, since they have no stars.

No generator output changes.

generation-output: unchanged
